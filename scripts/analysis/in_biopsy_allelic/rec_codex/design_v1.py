"""Complete-participant centering and fixed-schema Model A resampling.

Inputs must be synthetic or approved outcome-free designs until Gate B.
Raw z columns follow spec §3.3. Missing indicators follow the base columns
in their original order. Background and duplication are z columns 0 and 1.
"""
from dataclasses import dataclass, replace

import numpy as np

from likelihood_v1 import Model
from tests_v1 import ordinal_stage_probabilities, rng_for


@dataclass(frozen=True)
class Participant:
    individual: str
    cohort: str
    stage: int
    z: tuple
    e: float
    median_candidate_depth: float
    candidate_rows: int


class Design:
    def __init__(self, participants, rows, clusters=None, schema=None, gene_count=None):
        self.participants = tuple(sorted(participants, key=lambda p:p.individual))
        self.rows = tuple(rows)
        self.gene_count = gene_count if gene_count is not None else max((r.gene for r in rows),default=-1)+1
        self.by_id = {p.individual:p for p in self.participants}
        if len(self.by_id) != len(self.participants):
            raise ValueError('participant IDs must be unique')
        if any(r.individual not in self.by_id for r in rows):
            raise ValueError('row without participant')
        self.clusters = {p.individual:(clusters or {}).get(p.individual,p.individual) for p in self.participants}
        self.raw_z = np.asarray([p.z for p in self.participants],float)
        if self.raw_z.ndim != 2 or self.raw_z.shape[1] < 2:
            raise ValueError('raw z must start with background and duplication')
        self.cohort = np.asarray([p.cohort for p in self.participants])
        self.stage = np.asarray([p.stage for p in self.participants],float)
        self.stage_centered = self.stage.copy()
        for cohort in sorted(set(self.cohort)):
            mask = self.cohort == cohort
            self.stage_centered[mask] -= self.stage[mask].mean()
        missing = ~np.isfinite(self.raw_z)
        if schema is None:
            indicators = tuple(j for j in range(self.raw_z.shape[1]) if any(
                missing[self.cohort==c,j].any() and not missing[self.cohort==c,j].all()
                for c in set(self.cohort)))
            full = self._fill_center(indicators)
            keep = tuple(np.flatnonzero(np.any(full!=0,axis=0)).tolist())
            self.schema = {'indicators':indicators, 'keep':keep, 'raw_count':self.raw_z.shape[1]}
        else:
            self.schema = schema
            if self.raw_z.shape[1] != schema['raw_count']:
                raise ValueError('resample must preserve the raw z schema')
            full = self._fill_center(schema['indicators'])
        self.z = full[:,self.schema['keep']]
        self.background = full[:,0]
        self.duplication = full[:,1]

    def _fill_center(self, indicators):
        filled = self.raw_z.copy()
        missing = ~np.isfinite(filled)
        for cohort in sorted(set(self.cohort)):
            mask = self.cohort == cohort
            for j in range(filled.shape[1]):
                observed = mask & ~missing[:,j]
                filled[mask & missing[:,j],j] = filled[observed,j].mean() if observed.any() else 0.
        full = np.column_stack((filled,missing[:,indicators].astype(float)))
        for cohort in sorted(set(self.cohort)):
            mask = self.cohort == cohort
            constant=np.ptp(full[mask],axis=0)==0
            full[mask] -= full[mask].mean(axis=0)
            full[np.ix_(mask,constant)]=0.
        return full

    def model(self):
        index = {p.individual:j for j,p in enumerate(self.participants)}
        rows=[]
        for row in self.rows:
            j=index[row.individual]
            rows.append(replace(row,stage=float(self.stage_centered[j]),z=tuple(self.z[j]),
                                background=float(self.background[j]),duplication=float(self.duplication[j])))
        return Model(rows,self.gene_count)

    def with_stage(self, stage):
        stage=np.asarray(stage,int)
        if stage.shape != self.stage.shape or not np.isin(stage,(0,1,2)).all():
            raise ValueError('one stage in {0,1,2} per participant required')
        participants=[replace(p,stage=int(s)) for p,s in zip(self.participants,stage)]
        return Design(participants,self.rows,self.clusters,self.schema,self.gene_count)

    def stage_law(self):
        e=np.asarray([p.e for p in self.participants])
        depth=np.asarray([p.median_candidate_depth for p in self.participants])
        count=np.asarray([p.candidate_rows for p in self.participants])
        if np.any(~np.isfinite(e)) or np.any(e<=0) or np.any(count<0):
            raise ValueError('invalid error probability or candidate count')
        depth=np.where(count==0,0.,depth)
        if np.any(~np.isfinite(depth)) or np.any(depth<0):
            raise ValueError('invalid candidate median depth')
        x=np.column_stack((self.z,np.log(e),np.log1p(depth),np.log1p(count)))
        # Centering prevents a redundant common intercept with cohort thresholds.
        # Constant predictors carry no slope and are removed for this stage law.
        for cohort in set(self.cohort):
            mask=self.cohort==cohort
            constant=np.ptp(x[mask],axis=0)==0
            x[mask]-=x[mask].mean(axis=0)
            x[np.ix_(mask,constant)]=0.
        x=x[:,np.any(x!=0,axis=0)]
        return ordinal_stage_probabilities(self.stage,self.cohort,x)

    def stage_draws(self, draws, probabilities=None):
        if probabilities is None:
            fitted=self.stage_law()
            if not fitted['converged']:
                raise ValueError('stage law fit failed: '+fitted['message'])
            probabilities=fitted['probabilities']
        probabilities=np.asarray(probabilities,float)
        if probabilities.shape != (len(self.participants),3) or not np.allclose(probabilities.sum(axis=1),1):
            raise ValueError('invalid stage-law probabilities')
        for b in range(draws):
            u=rng_for(1,b).random(len(self.participants))
            yield np.sum(u[:,None]>np.cumsum(probabilities,axis=1),axis=1)

    def cohort_clusters(self):
        groups={}
        home={}
        for p in self.participants:
            cluster=self.clusters[p.individual]
            # Participants are sorted: the first member binds the v1.2 stratum.
            home.setdefault(cluster,p.cohort)
            groups.setdefault(home[cluster],{}).setdefault(cluster,[]).append(p.individual)
        return {c:[v[k] for k in sorted(v)] for c,v in sorted(groups.items())}

    def cluster_indices(self, b):
        rng=rng_for(2,b)
        return {c:rng.integers(0,len(v),len(v)) for c,v in self.cohort_clusters().items()}

    def resample(self, indices):
        groups=self.cohort_clusters()
        if set(indices) != set(groups):
            raise ValueError('resample must supply every cohort')
        participants=[];rows=[];clusters={}
        rows_by_id={p.individual:[] for p in self.participants}
        for row in self.rows:
            rows_by_id[row.individual].append(row)
        for cohort in sorted(groups):
            picks=np.asarray(indices[cohort],int)
            if len(picks) != len(groups[cohort]):
                raise ValueError('resample must keep the cohort cluster count')
            for copy, pick in enumerate(picks):
                for original in groups[cohort][pick]:
                    individual=f'{original}__copy{copy:06d}'
                    clusters[individual]=f'{self.clusters[original]}__copy{copy:06d}'
                    participants.append(replace(self.by_id[original],individual=individual))
                    rows.extend(replace(row,individual=individual) for row in rows_by_id[original])
        return Design(participants,rows,clusters,self.schema,self.gene_count)

    def stage_r_squared(self,adjusted=True):
        result={}
        for cohort in sorted(set(self.cohort)):
            mask=self.cohort==cohort
            y=self.stage_centered[mask];x=self.z[mask]
            total=float(y@y)
            if total<=0:
                result[cohort]=float('nan')
                continue
            beta,_,rank,_=np.linalg.lstsq(x,y,rcond=1e-12)
            residual=y-x@beta
            raw=float(np.clip(1-residual@residual/total,0,1))
            n=int(mask.sum())
            result[cohort]=(float(max(0.,1-(1-raw)*(n-1)/(n-rank-1)))
                            if adjusted and n>rank+1 else
                            float('nan') if adjusted else raw)
        return result
