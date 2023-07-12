#conda install pytorch==1.7.1 torchvision==0.8.2 torchaudio==0.7.2 -c pytorch

~/miniconda3/bin/python
## 1) login through gpu node
srun --pty  --cpus-per-task=1 --job-name=interact --ntasks=1 --nodes=1 --partition=gpu_p --time=12:00:00 --mem=2GB /bin/bash -l
srun --pty  --cpus-per-task=1 --job-name=interact --ntasks=1 --nodes=1 --partition=gpu_p --time=1:00:00 --gres=gpu:P100:1 --mem=13G /bin/bash -l
srun --pty  --cpus-per-task=1 --job-name=interact --ntasks=1 --nodes=1 --partition=gpu_p --time=1:00:00 --gres=gpu:A100:1 --mem=13G /bin/bash -l

## 2) conda activate
conda activate pytorch

## 3) module load
module load CUDA/11.1.1-GCC-10.2.0

cd /scratch/sb14489/8.ML_ACR/DeepFormer_Ex/DeepFormer/maize_code


from selene_sdk.utils import load_path
from selene_sdk.utils import parse_configs_and_run

configs = load_path("./DanQ_YAML.yml")

parse_configs_and_run(configs,lr=0.0005)
