import importlib,importlib.metadata,json,pathlib,subprocess,inspect
import torch
import numpy as np
root=pathlib.Path('Analysis/ATAC/Regulatory_Methods/envs/installation-20261023')
modules=['LingerGRN.preprocess','LingerGRN.LINGER_tr','LingerGRN.LL_net','LingerGRN.pseudo_bulk','scanpy','shap','pybedtools']
results={}
for name in modules:
 try:
  obj=importlib.import_module(name);results[name]={'status':'passed','source':obj.__file__}
 except Exception as exc:results[name]={'status':'failed','exception':repr(exc)}
torch.set_num_threads(1);torch.manual_seed(42)
x=torch.tensor([[1.,2.],[3.,4.]])
model=torch.nn.Linear(2,1);loss=model(x).square().mean();loss.backward()
results['torch']={'version':torch.__version__,'cuda_available':torch.cuda.is_available(),'compiled_cuda':torch.version.cuda,'forward_backward_finite':bool(torch.isfinite(loss))}
results['bedtools']=subprocess.check_output(['bedtools','--version'],text=True).strip()
(root/'linger-smoke.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
if any(v.get('status')=='failed' for v in results.values() if isinstance(v,dict)):raise SystemExit(1)
