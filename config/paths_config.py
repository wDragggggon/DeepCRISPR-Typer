import os
from pathlib import Path

# 项目根目录
PROJECT_ROOT = Path(__file__).parent.parent

# 输入输出目录
DATA_DIR = PROJECT_ROOT / "data"
INPUT_DIR = DATA_DIR / "input"
OUTPUT_DIR = DATA_DIR / "output"

# 模型目录
MODELS_DIR = DATA_DIR / "models"

# 外部工具路径 - 用户需要根据自己的系统配置这些路径
# MinCED: https://github.com/ctSkennerton/minced
# Prodigal: http://prodigal.ornl.gov/
# HMMER: http://hmmer.org/
MINCED_PATH = "minced.jar"  # 请替换为您的MinCED jar文件实际路径
PRODIGAL_PATH = "prodigal"  # 请替换为您的Prodigal可执行文件路径（通常在PATH中）
HMMSCAN_PATH = "hmmscan"    # 请替换为您的HMMER hmmscan可执行文件路径（通常在PATH中）

# 现有模型路径 - 用户需要下载预训练模型并放置在相应位置
# CNN-Att模型来自: CRISPRclassify-CNN-Att (请参考相关论文获取)
CNN_ATT_MODEL_PATH = str(MODELS_DIR / "cnn_att_large.pth")  # 请将预训练模型放在此路径

# TEMC-Cas模型来自: TEMC-Cas (请参考相关论文获取)
TEMC_CAS_MODEL_PATH = str(MODELS_DIR / "temc_cas_best_model.pth")  # 请将预训练模型放在此路径

# CRISPRCasTyper数据路径 - 用户需要下载相关数据文件
# 这些文件通常来自CRISPRCasTyper项目: https://github.com/Russel88/CRISPRCasTyper
HMM_DB_PATH = str(DATA_DIR / "Profiles")  # 请将HMM数据库放在此目录
ADAPTATION_RULES_PATH = str(DATA_DIR / "adaptation.json")  # 请将适应规则文件放在此路径
INTERFERENCE_RULES_PATH = str(DATA_DIR / "interference.json")  # 请将干扰规则文件放在此路径

# 专家规则权重矩阵路径
CAS_SCORING_PATH = PROJECT_ROOT / "CasScoring.csv"

# 创建必要的目录
for dir_path in [INPUT_DIR, OUTPUT_DIR, MODELS_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)

# 验证关键路径是否存在
def validate_paths():
    """验证所有必需的路径是否存在"""
    required_paths = {
        "MinCED": MINCED_PATH,
        "Prodigal": PRODIGAL_PATH,
        "HMMER": HMMSCAN_PATH,
        "CNN-Att Model": CNN_ATT_MODEL_PATH,
        "HMM Database": HMM_DB_PATH,
        "Adaptation Rules": ADAPTATION_RULES_PATH,
        "Interference Rules": INTERFERENCE_RULES_PATH,
        "CasScoring Matrix": CAS_SCORING_PATH
    }
    
    # TEMC-Cas模型是可选的
    if os.path.exists(TEMC_CAS_MODEL_PATH):
        required_paths["TEMC-Cas Model"] = TEMC_CAS_MODEL_PATH
    else:
        print("⚠️  TEMC-Cas Model not found at expected location")
    
    missing_paths = []
    for name, path in required_paths.items():
        if not os.path.exists(path):
            missing_paths.append(f"{name}: {path}")
    
    if missing_paths:
        print("⚠️  Missing required paths:")
        for path in missing_paths:
            print(f"  - {path}")
        return False
    else:
        print("✅ All required paths are valid")
        return True