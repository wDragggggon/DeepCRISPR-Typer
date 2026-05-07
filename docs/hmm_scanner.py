import os
import subprocess
import tempfile
import re
from typing import List, Dict, Optional, Tuple
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed # 在文件顶部导入
import pandas as pd

from config.paths_config import HMM_DB_PATH, HMMSCAN_PATH, CAS_SCORING_PATH
from config.model_config import HMM_CONFIG
from src.utils import run_command, load_fasta_sequences


# 【新增】：必须放在类外面，才能被多进程序列化
def _hmm_scan_worker(task_data):
    """
    单个蛋白扫描的独立进程任务
    task_data: 包含所有必要参数的字典
    """
    import os
    import subprocess
    from src.utils import run_command # 确保 worker 能找到这个函数
    
    # 解包参数
    protein_id = task_data['protein_id']
    sequence = task_data['sequence']
    relevant_hmms = task_data['relevant_hmms']
    hmm_file_index = task_data['hmm_file_index']
    hmmscan_path = task_data['hmmscan_path']
    evalue_threshold = task_data['evalue_threshold']
    bitscore_threshold = task_data['bitscore_threshold']
    domain_evalue_threshold = task_data['domain_evalue_threshold']

    # 注意：这里需要重新实现一些辅助逻辑，因为 Worker 无法访问 self
    # 1. 创建临时 HMM DB (由于逻辑复杂，建议直接在这里快速实现简易版)
    import tempfile
    try:
        # 创建 HMM
        with tempfile.NamedTemporaryFile(mode='w', suffix='.hmm', delete=False) as tmp_hmm:
            for hmm_name in relevant_hmms:
                if hmm_name in hmm_file_index:
                    with open(hmm_file_index[hmm_name], 'r') as f:
                        tmp_hmm.write(f.read())
                        tmp_hmm.write("\n")
            tmp_hmm_name = tmp_hmm.name
        
        # 执行 hmmpress
        subprocess.run(["hmmpress", "-f", tmp_hmm_name], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        # 2. 创建临时 FASTA
        with tempfile.NamedTemporaryFile(mode='w', suffix='.faa', delete=False) as tmp_fa:
            tmp_fa.write(f">{protein_id}\n{sequence}\n")
            tmp_fa_name = tmp_fa.name

        # 3. 运行 hmmscan
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp_out:
            output_file = tmp_out.name

        cmd = [hmmscan_path, "--tblout", output_file, "--cpu", "1", tmp_hmm_name, tmp_fa_name]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # 4. 解析结果 (复用简单的解析逻辑)
        results = []
        if os.path.exists(output_file):
            with open(output_file, 'r') as f:
                for line in f:
                    if line.startswith('#'): continue
                    fields = line.strip().split()
                    if len(fields) < 18: continue
                    e_val, score = float(fields[4]), float(fields[5])
                    dom_e_val, dom_score = float(fields[12]), float(fields[13])
                    
                    if e_val <= evalue_threshold and score >= bitscore_threshold:
                        results.append({
                            'protein_id': fields[2],
                            'hmm_family': fields[0],
                            'score': score,
                            'e_value': e_val
                        })
        
        # 清理
        for f in [tmp_hmm_name, tmp_fa_name, output_file] + [tmp_hmm_name + e for e in ['.h3m','.h3i','.h3f','.h3p']]:
            if os.path.exists(f): os.unlink(f)
            
        return results
    except Exception as e:
        return []

class HMMScanner:
    """HMM扫描器，用于Cas蛋白家族匹配 - 实现定向比对策略"""
    
    def __init__(self, hmm_db_path: str = None, hmmscan_path: str = None):
        self.hmm_db_path = hmm_db_path or HMM_DB_PATH
        self.hmmscan_path = hmmscan_path or HMMSCAN_PATH
        self.evalue_threshold = HMM_CONFIG['evalue_threshold']
        self.bitscore_threshold = HMM_CONFIG['bitscore_threshold']
        self.domain_evalue_threshold = HMM_CONFIG['domain_evalue_threshold']
        
        # 加载专家规则权重矩阵
        self.cas_scoring_matrix = self._load_cas_scoring_matrix(CAS_SCORING_PATH)
        
        # 构建HMM文件索引和家族映射
        self.hmm_file_index = self._build_hmm_index()
        self.family_to_hmms = self._build_family_to_hmms_mapping()
    
    def _load_cas_scoring_matrix(self, scoring_path: str) -> pd.DataFrame:
        """加载专家规则权重矩阵"""
        try:
            if not Path(scoring_path).exists():
                print(f"Warning: CasScoring matrix not found at {scoring_path}")
                return None
            
            df = pd.read_csv(scoring_path, index_col=0, na_values=['', ' ', 'NaN'])
            df = df.fillna(0.0)
            return df
        except Exception as e:
            print(f"Failed to load CasScoring matrix: {e}")
            return None
    
    def _build_hmm_index(self) -> Dict[str, str]:
        """构建HMM文件名到完整路径的索引"""
        hmm_index = {}
        if os.path.exists(self.hmm_db_path):
            hmm_files = []
            for root, dirs, files in os.walk(self.hmm_db_path):
                for file in files:
                    if file.endswith('.hmm'):
                        full_path = os.path.join(root, file)
                        hmm_name = os.path.splitext(file)[0]  # 移除.hmm扩展名
                        hmm_index[hmm_name] = full_path
                        hmm_files.append(hmm_name)
            print(f"✅ Indexed {len(hmm_files)} HMM models")
        else:
            print(f"Warning: HMM database path not found: {self.hmm_db_path}")
        return hmm_index
    
    def _build_family_to_hmms_mapping(self) -> Dict[str, List[str]]:
        """构建Cas蛋白家族到相关HMM模型的映射"""
        family_mapping = {}
        
        # 定义Cas家族关键词到亚型的映射规则
        family_keywords = {
            'Cas1': ['I', 'II', 'III', 'IV', 'V', 'VI'],
            'Cas2': ['I', 'II', 'III', 'IV', 'V', 'VI'],
            'Cas3': ['I'],
            'Cas4': ['I', 'II', 'V'],
            'Cas5': ['I'],
            'Cas6': ['I', 'III', 'IV'],
            'Cas7': ['I'],
            'Cas8': ['I'],
            'Cas9': ['II'],
            'Cas10': ['III'],
            'Cas12': ['V'],
            'Cas13': ['VI'],
            'Cmr': ['III'],  # Cmr1, Cmr3, Cmr4, etc.
            'Csm': ['III'],  # Csm2, Csm3, Csm4, etc.
            'Csf': ['IV'],   # Csf1, Csf2, Csf3, Csf4
            'Csn': ['II'],   # Csn2
        }
        
        # 为每个HMM文件确定其相关的Cas家族
        for hmm_name in self.hmm_file_index.keys():
            # 提取HMM名称中的家族信息
            detected_families = set()
            
            # 检查是否包含特定家族关键词
            for family, keywords in family_keywords.items():
                if family in hmm_name:
                    detected_families.add(family)
            
            # 特殊处理：如果包含数字后跟下划线，可能是Cas+数字的模式
            if re.search(r'Cas\d+', hmm_name):
                # 提取Cas数字
                cas_match = re.search(r'Cas(\d+)', hmm_name)
                if cas_match:
                    cas_num = cas_match.group(1)
                    if cas_num == '1':
                        detected_families.add('Cas1')
                    elif cas_num == '2':
                        detected_families.add('Cas2')
                    elif cas_num == '3':
                        detected_families.add('Cas3')
                    elif cas_num == '9':
                        detected_families.add('Cas9')
                    elif cas_num == '10':
                        detected_families.add('Cas10')
                    elif cas_num == '12':
                        detected_families.add('Cas12')
                    elif cas_num == '13':
                        detected_families.add('Cas13')
            
            # 将HMM添加到对应的家族列表中
            for family in detected_families:
                if family not in family_mapping:
                    family_mapping[family] = []
                family_mapping[family].append(hmm_name)
        
        print(f"✅ Built family-to-HMM mapping for {len(family_mapping)} families")
        return family_mapping
    
    # def scan_proteins_with_families(self, protein_fasta_file: str, cas_families: List[Dict]) -> List[Dict]:
    #     """对蛋白序列进行定向HMM扫描"""
    #     # ... (保留函数前面的验证逻辑和循环开始部分) ...
    #     if not os.path.exists(protein_fasta_file): return []
    #     if not self.hmm_file_index: return []
        
    #     protein_sequences = load_fasta_sequences(protein_fasta_file)
    #     all_results = []
        
    #     for cas_pred in cas_families:
    #         protein_id = cas_pred['protein_id']
    #         predicted_family = cas_pred['predicted_family']
            
    #         if protein_id not in protein_sequences: continue
            
    #         relevant_hmms = self._get_relevant_hmms_for_family(predicted_family)
    #         if not relevant_hmms: relevant_hmms = self._fuzzy_match_hmms(predicted_family)
            
    #         if relevant_hmms:
    #             temp_hmm_db = self._create_temp_hmm_db(relevant_hmms)
    #             if temp_hmm_db:
    #                 try:
    #                     temp_fasta = self._create_temp_fasta_for_protein(protein_id, protein_sequences[protein_id])
    #                     if temp_fasta:
    #                         try:
    #                             results = self._scan_single_protein(temp_fasta, temp_hmm_db, protein_id)
    #                             all_results.extend(results)
    #                         finally:
    #                             if os.path.exists(temp_fasta): os.unlink(temp_fasta)
    #                 finally:
    #                     # 【关键修复】：不仅删除 .hmm，还要删除 hmmpress 生成的 4 个二进制索引文件
    #                     if os.path.exists(temp_hmm_db):
    #                         os.unlink(temp_hmm_db)
    #                     for ext in ['.h3m', '.h3i', '.h3f', '.h3p']:
    #                         if os.path.exists(temp_hmm_db + ext):
    #                             os.unlink(temp_hmm_db + ext)
    #     return all_results
    
    def scan_proteins_with_families(self, protein_fasta_file: str, cas_families: List[Dict]) -> List[Dict]:
        if not os.path.exists(protein_fasta_file) or not self.hmm_file_index:
            return []
        
        protein_sequences = load_fasta_sequences(protein_fasta_file)
        tasks = []
        
        # 准备任务数据包
        for cas_pred in cas_families:
            prot_id = cas_pred['protein_id']
            if prot_id not in protein_sequences: continue
            
            relevant_hmms = self._get_relevant_hmms_for_family(cas_pred['predicted_family'])
            if not relevant_hmms: relevant_hmms = self._fuzzy_match_hmms(cas_pred['predicted_family'])
            
            if relevant_hmms:
                tasks.append({
                    'protein_id': prot_id,
                    'sequence': protein_sequences[prot_id],
                    'relevant_hmms': list(set(relevant_hmms)),
                    'hmm_file_index': self.hmm_file_index,
                    'hmmscan_path': self.hmmscan_path,
                    'evalue_threshold': self.evalue_threshold,
                    'bitscore_threshold': self.bitscore_threshold,
                    'domain_evalue_threshold': self.domain_evalue_threshold
                })

        all_results = []
        # 使用进程池火力全开
        with ProcessPoolExecutor(max_workers=32) as executor:
            future_to_task = {executor.submit(_hmm_scan_worker, t): t for t in tasks}
            for future in as_completed(future_to_task):
                res = future.result()
                if res: all_results.extend(res)
                    
        return all_results
    
    def _get_relevant_hmms_for_family(self, family: str) -> List[str]:
        """获取指定家族相关的HMM模型列表"""
        # 直接匹配
        if family in self.family_to_hmms:
            return self.family_to_hmms[family]
        
        # 模糊匹配：检查family是否是某个已知家族的子字符串
        for known_family in self.family_to_hmms:
            if family.startswith(known_family) or known_family.startswith(family):
                return self.family_to_hmms[known_family]
        
        return []
    
    def _fuzzy_match_hmms(self, family: str) -> List[str]:
        """模糊匹配HMM模型"""
        matched_hmms = []
        family_lower = family.lower()
        
        for hmm_name in self.hmm_file_index.keys():
            if family_lower in hmm_name.lower():
                matched_hmms.append(hmm_name)
        
        return matched_hmms
    
    def _create_temp_hmm_db(self, hmm_names: List[str]) -> str:
        """创建包含指定HMM模型的临时数据库文件，并执行hmmpress建立索引"""
        if not hmm_names:
            return None
            
        # 列表去重
        unique_hmm_names = list(set(hmm_names))
        
        try:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.hmm', delete=False) as tmp:
                tmp_name = tmp.name
                
            with open(tmp_name, 'w') as f_out:
                for hmm_name in unique_hmm_names:
                    if hmm_name in self.hmm_file_index:
                        hmm_path = self.hmm_file_index[hmm_name]
                        process = subprocess.run(
                            ["hmmconvert", hmm_path], 
                            capture_output=True, 
                            text=True
                        )
                        if process.returncode != 0:
                            print(f"⚠️ hmmconvert 警告 ({hmm_name}): {process.stderr}")
                        else:
                            # 【终极必杀技】：强制重写 HMM 内部的 NAME 属性，将其替换为唯一的文件名！
                            hmm_text = process.stdout
                            hmm_text = re.sub(r'^NAME\s+.*$', f'NAME  {hmm_name}', hmm_text, flags=re.MULTILINE)
                            f_out.write(hmm_text)
                            f_out.write("\n")
            
            # 执行 hmmpress
            process = subprocess.run(
                ["hmmpress", "-f", tmp_name], 
                capture_output=True, text=True
            )
            
            if process.returncode != 0:
                print(f"⚠️ hmmpress failed! Error: {process.stderr}")
                os.unlink(tmp_name)
                return None
                
            return tmp_name
            
        except Exception as e:
            print(f"Failed to create temporary HMM database: {e}")
            if 'tmp_name' in locals() and os.path.exists(tmp_name):
                os.unlink(tmp_name)
            return None
    
    def _create_temp_fasta_for_protein(self, protein_id: str, sequence: str) -> str:
        """为单个蛋白序列创建临时FASTA文件"""
        try:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.faa', delete=False) as tmp:
                tmp.write(f">{protein_id}\n{sequence}\n")
                return tmp.name
        except Exception as e:
            print(f"Failed to create temporary FASTA file: {e}")
            return None
    
    def _scan_single_protein(self, protein_fasta: str, hmm_db: str, protein_id: str) -> List[Dict]:
        """对单个蛋白序列进行HMM扫描"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp_out:
            output_file = tmp_out.name
        
        try:
            cmd = [
                self.hmmscan_path,
                "--tblout", output_file,
                "--cpu", "1",
                # 【核心修复 2】：删除了 "--cut_ga" 参数，交由我们自己的 Python 代码进行阈值控制！
                hmm_db,
                protein_fasta
            ]
            
            stdout, stderr, returncode = run_command(cmd)
            
            if returncode != 0:
                print(f"HMM scan failed for {protein_id}: {stderr}")
                return []
            
            results = self._parse_hmmscan_output(output_file, protein_id)
            return results
            
        finally:
            if os.path.exists(output_file):
                os.unlink(output_file)
    
    def _parse_hmmscan_output(self, tblout_file: str, expected_protein_id: str = None) -> List[Dict]:
        """解析hmmscan的tblout格式输出"""
        results = []
        
        if not os.path.exists(tblout_file):
            return results
        
        with open(tblout_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split()
                if len(fields) < 18:
                    continue
                
                target_name = fields[0]      # HMM名称
                query_name = fields[2]       # 查询序列名称
                e_value = float(fields[4])   # 全序列E值
                score = float(fields[5])     # 全序列得分
                domain_e_value = float(fields[12])  # 域E值
                domain_score = float(fields[13])    # 域得分
                
                # 应用阈值过滤
                if (e_value <= self.evalue_threshold and 
                    score >= self.bitscore_threshold and
                    domain_e_value <= self.domain_evalue_threshold):
                    
                    # 提取亚型信息
                    subtype = self._extract_subtype_from_hmm(target_name)
                    
                    results.append({
                        'protein_id': query_name,
                        'hmm_family': target_name,
                        'subtype': subtype,
                        'e_value': e_value,
                        'score': score,
                        'domain_e_value': domain_e_value,
                        'domain_score': domain_score
                    })
        
        return results
    
    def _extract_subtype_from_hmm(self, hmm_name: str) -> str:
        """从HMM名称中提取亚型信息"""
        # 尝试从HMM名称中提取亚型
        # 例如: Cas12a_0_CAS-V-A.hmm -> V-A
        #       Csf4_2_IVA1.hmm -> IV-A1
        
        # 查找CAS-前缀的模式
        cas_match = re.search(r'CAS-([IVX\-A-Z0-9]+)', hmm_name)
        if cas_match:
            subtype = cas_match.group(1).replace('-', '')
            # 标准化亚型格式
            if subtype.startswith('I') and not subtype.startswith('II') and not subtype.startswith('III') and not subtype.startswith('IV') and not subtype.startswith('V') and not subtype.startswith('VI'):
                # Type I 的子类型
                if len(subtype) > 1:
                    return f"I-{subtype[1:]}"
                else:
                    return "I"
            elif subtype.startswith('II'):
                if len(subtype) > 2:
                    return f"II-{subtype[2:]}"
                else:
                    return "II"
            elif subtype.startswith('III'):
                if len(subtype) > 3:
                    return f"III-{subtype[3:]}"
                else:
                    return "III"
            elif subtype.startswith('IV'):
                if len(subtype) > 2:
                    return f"IV-{subtype[2:]}"
                else:
                    return "IV"
            elif subtype.startswith('V'):
                if len(subtype) > 1:
                    return f"V-{subtype[1:]}"
                else:
                    return "V"
            elif subtype.startswith('VI'):
                if len(subtype) > 2:
                    return f"VI-{subtype[2:]}"
                else:
                    return "VI"
        
        # 查找罗马数字模式
        roman_match = re.search(r'([IVX]+)([A-Z0-9]*)', hmm_name)
        if roman_match:
            roman_part = roman_match.group(1)
            subtype_part = roman_match.group(2)
            if subtype_part:
                return f"{roman_part}-{subtype_part}"
            else:
                return roman_part
        
        return hmm_name  # 如果无法解析，返回原始名称
    
    # 【关键修改】：将函数名和逻辑改为按蛋白生成向量
    def calculate_protein_propensity_vectors(self, hmm_results: List[Dict]) -> Dict[str, Dict[str, float]]:
        """计算每个候选蛋白 Sj 的亚型倾向性向量 SP,j"""
        protein_vectors = {}
        
        # 获取所有可能的亚型列表作为向量的维度
        subtypes = list(self.cas_scoring_matrix.columns) if self.cas_scoring_matrix is not None else []
        
        for result in hmm_results:
            prot_id = result['protein_id']
            hmm_family = result['hmm_family']
            bit_score = result['score'] # 提取 Bit Score 阈值过滤后的得分
            
            # 初始化该蛋白的零向量
            if prot_id not in protein_vectors:
                if self.cas_scoring_matrix is not None:
                    protein_vectors[prot_id] = {col: 0.0 for col in subtypes}
                else:
                    protein_vectors[prot_id] = {}

            # 向量累加：Weight * BitScore
            if self.cas_scoring_matrix is not None and hmm_family in self.cas_scoring_matrix.index:
                weights = self.cas_scoring_matrix.loc[hmm_family]
                for subtype, weight in weights.items():
                    if pd.notna(weight) and weight != 0:
                        protein_vectors[prot_id][subtype] += float(weight) * bit_score
            else:
                # Fallback：如果没有权重矩阵，直接累加 Bit Score
                subtype = result.get('subtype', hmm_family)
                protein_vectors[prot_id][subtype] = protein_vectors[prot_id].get(subtype, 0.0) + bit_score
                
        return protein_vectors