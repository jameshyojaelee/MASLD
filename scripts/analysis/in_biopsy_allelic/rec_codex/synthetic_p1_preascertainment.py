"""Synthetic Model A null generator; no fit, resampling test, or real input.

The two regimes share pre-read design, genotypes, and count uniforms. All
candidate tags receive reads before depth-only tag selection and the integer
het rule. Within a gene the two tags represent the same synthetic haplotypes,
with reversed allele coding; this is not a model of real linkage or read overlap.
"""
from dataclasses import dataclass, replace

import numpy as np
from scipy.special import expit, logit
from scipy.stats import betabinom

from design_v1 import Design, Participant
from likelihood_v1 import EPS, Row


SEED = 20260930
ALPHA = np.array([.7, -.6, 1.1, -.3, .5, -.9, .8, -.4])
RAW_Z_NAMES = ('background_log_error', 'duplication', 'ffpe', 'sex_female',
               'age10', 'axis_1', 'axis_2')
RULES = ((20, 3, 15), (10, 2, 10), (10, 1, 5))
REGIMES = {'nuisance_off': (0., 0.), 'nuisance_on': (.02, .10)}


@dataclass(frozen=True)
class SyntheticNull:
    design: Design
    candidate_tags: tuple
    chosen_candidates: tuple
    selected_keys: tuple
    preselection_stage_probabilities: np.ndarray
    generating_parameters: np.ndarray
    metadata: dict


def generate(regime, *, seed=SEED, replicate=0, participants_per_cohort=60):
    """Return a complete synthetic participant design and selected RNA rows.

    Stage is generated with a correctly specified cohort-threshold ordinal
    logit using age and sex. Read generation uses complete-participant realized
    cohort centering, exactly as Design does. Consequently the preselection
    stage probabilities are NOT an oracle law conditional on the inclusion
    pattern; cohort centering couples that pattern across participants.
    """
    if regime not in REGIMES:
        raise ValueError('Expected nuisance_off or nuisance_on.')
    if int(participants_per_cohort) != participants_per_cohort or participants_per_cohort < 12:
        raise ValueError('At least twelve integer participants per cohort required.')
    if int(replicate) != replicate or replicate < 0:
        raise ValueError('Nonnegative integer replicate required.')
    # Independent streams avoid nuisance-dependent RNG consumption changing design.
    streams = [np.random.default_rng(np.random.SeedSequence([seed, replicate, j]))
               for j in range(4)]
    rd, rs, rg, ra = streams
    n = 3 * participants_per_cohort
    cohort = np.repeat(np.arange(3), participants_per_cohort)
    error = np.clip(np.exp(rd.normal(np.log(.002), .3, n)), 1e-4, .02)
    duplication = rd.uniform(.05, .4, n)
    sex = rd.integers(0, 2, n).astype(float)
    age = rd.normal(5.2, .8, n)
    axes = rd.normal(0., .5, (n, 2))
    raw_z = np.column_stack((np.log(error), duplication, np.zeros(n), sex, age, axes))
    centered = raw_z.copy()
    for c in range(3):
        centered[cohort == c] -= centered[cohort == c].mean(axis=0)
    predictor = .9 * centered[:, 4] + .3 * centered[:, 3]
    thresholds = np.array([[-.8, .6], [-.6, .8], [-.7, .7]])
    cumulative = expit(thresholds[cohort] - predictor[:, None])
    probabilities = np.column_stack((cumulative[:, 0], cumulative[:, 1]-cumulative[:, 0],
                                     1-cumulative[:, 1]))
    stage = np.sum(rs.random(n)[:, None] > np.cumsum(probabilities, axis=1), axis=1)
    # Depth depends on recorded age; random dropout is independent of allele reads.
    depths = np.clip(np.rint(np.exp(rd.normal(np.log(35), .7, (n, 8, 2))
                           + .12 * centered[:, 4, None, None])), 0, 500).astype(int)
    depths[rd.random((n, 8)) < .08] = 0
    # A metadata-fixed negative-control group has no candidates and remains in F.
    for c in range(3):
        depths[np.flatnonzero(cohort == c)[:2]] = 0
    omega = (.24, .095, .08)
    people, templates, chosen = [], [], []
    for i in range(n):
        c = int(cohort[i]); label = f'C{c+1}'
        identity = f'{label}_P{i:04d}'
        rule = RULES[c]
        candidate_depths = []
        for g in range(8):
            tags = []
            for t in range(2):
                orientation = 1 if (g+t) % 2 == 0 else -1
                # Reverse genomic allele coding, retaining one shared oriented genotype.
                ref_frequency, alt_frequency = ((.7, .3) if orientation == 1 else (.3, .7))
                row = Row(g, identity, label, int(depths[i, g, t]), 0,
                          ref_frequency, alt_frequency, 0., float(error[i]), .0439,
                          omega[c], orientation, 0., tuple(raw_z[i]), 0., 0., *rule)
                tags.append(row); templates.append((i, g, t, row))
            t = int(np.argmax(depths[i, g]))  # ties to the lower synthetic tag position
            chosen.append((i, g, t, tags[t]))
            if tags[t].n >= rule[0]:
                candidate_depths.append(tags[t].n)
        people.append(Participant(identity, label, int(stage[i]), tuple(raw_z[i]), float(error[i]),
                                  float(np.median(candidate_depths)) if candidate_depths else 0.,
                                  len(candidate_depths)))
    pre_design = Design(people, [r for _, _, _, r in chosen], gene_count=8)
    model = pre_design.model()
    p = np.zeros(model.size)
    p[:8] = ALPHA
    p[8:16] = logit((.02-EPS)/(1-2*EPS))
    # Nonzero recorded-covariate dependence is present in both null regimes.
    p[17:17+model.Z] = np.array([.02, -.10, 0., .03, .01, .02, -.02])[list(pre_design.schema['keep'])]
    omega_s, rho_s = REGIMES[regime]
    p[17+model.Z] = omega_s
    p[-3:] = (rho_s, .05, -.05)
    genotypes = rg.choice(3, size=(n, 8), p=(.42, .49, .09))
    uniforms = ra.random((n, 8, 2))
    generated = {}
    participant_index = {person.individual:j for j, person in enumerate(pre_design.participants)}
    class_selected = np.zeros(3, dtype=int)
    class_candidates = np.bincount(genotypes.ravel(), minlength=3)
    for i, g, t, row in templates:
        j = participant_index[row.individual]
        current = replace(row, stage=float(pre_design.stage_centered[j]),
                          z=tuple(pre_design.z[j]),
                          background=float(pre_design.background[j]),
                          duplication=float(pre_design.duplication[j]))
        eta, r, _, _ = model.predictors(current, p, need_hessian=False)
        cls = int(genotypes[i, g])
        if cls == 0:
            mean = expit(eta-current.orientation*current.omega)
            rho = EPS+(1-2*EPS)*expit(r)
        else:
            mean = current.e if cls == 1 else 1-current.e
            rho = current.phi
        concentration = (1-rho)/rho
        count = int(betabinom.ppf(uniforms[i, g, t], current.n,
                                 mean*concentration, (1-mean)*concentration))
        if not 0 <= count <= current.n:
            raise ValueError('Generated count is outside its RNA read depth.')
        generated[i, g, t] = replace(row, a=count)
    included, keys = [], []
    chosen_read_rows = []
    for i, g, t, row in chosen:
        read_row = generated[i, g, t]
        chosen_read_rows.append(read_row)
        if read_row.a in read_row.allowed:
            included.append(read_row); keys.append((row.individual, g, t))
            class_selected[genotypes[i, g]] += 1
    final = Design(people, included, schema=pre_design.schema, gene_count=8)
    metadata = {'seed':seed, 'replicate':replicate, 'regime':regime,
                'participants':n, 'genes':8, 'tags_per_gene':2, 'raw_z_names':RAW_Z_NAMES,
                'kappa':0., 'tau':0., 'omega_S':omega_s, 'rho_S':rho_s,
                'stage_law':'cohort ordered thresholds; common age and sex slopes .9 and .3',
                'stage_counts':{f'C{c+1}':np.bincount(stage[cohort == c], minlength=3).tolist()
                                for c in range(3)},
                'candidate_genotype_class_counts_H_R_A':class_candidates.tolist(),
                'selected_genotype_class_counts_H_R_A':class_selected.tolist(),
                'selected_rows':len(included),
                'zero_selected_participants':sum(not any(r.individual == person.individual for r in included)
                                                 for person in people),
                'zero_candidate_participants':sum(person.candidate_rows == 0 for person in people),
                'preselection_probabilities_are_selection_conditional':False,
                'fits_or_tests_performed':False,
                'limits':'Synthetic linked tag coding; independent read draws across tags/genes; '
                         'fixed candidate depth law; no real RNA design or read overlap. '
                         'No P1 size, power, or oracle conditional stage-law claim.'}
    return SyntheticNull(final, tuple(generated[i,g,t] for i,g,t,_ in templates),
                         tuple(chosen_read_rows), tuple(keys), probabilities, p, metadata)
