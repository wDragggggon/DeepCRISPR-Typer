import json
import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Optional
from pathlib import Path

from config.paths_config import ADAPTATION_RULES_PATH, INTERFERENCE_RULES_PATH, CAS_SCORING_PATH
from config.model_config import RULE_CONFIG

class RuleBasedClassifier:
    """基于规则的CRISPR-Cas系统分类器"""
    
    def __init__(self, adaptation_rules_path: str = None, interference_rules_path: str = None, cas_scoring_path: str = None):
        self.adaptation_rules_path = adaptation_rules_path or ADAPTATION_RULES_PATH
        self.interference_rules_path = interference_rules_path or INTERFERENCE_RULES_PATH
        self.cas_scoring_path = cas_scoring_path or CAS_SCORING_PATH
        
        self.adaptation_rules = self._load_rules(self.adaptation_rules_path)
        self.interference_rules = self._load_rules(self.interference_rules_path)
        self.cas_scoring_matrix = self._load_cas_scoring_matrix(self.cas_scoring_path)
        
        self.confidence_threshold = RULE_CONFIG['final_confidence_threshold']
        
    def _aggregate_locus_score(self, protein_vectors: Dict[str, Dict[str, float]]) -> Tuple[Dict, str]:
        """执行论文逻辑：计算基因座级别向量 S_Locus = \sum S_{P,j}"""
        if not protein_vectors:
            return {'adaptation_score': 0.0, 'interference_score': 0.0}, ""
            
        # 提取所有亚型维度
        subtypes = list(self.cas_scoring_matrix.columns)
        s_locus = {subtype: 0.0 for subtype in subtypes}
        
        # 将基因座内所有蛋白的向量 S_{P,j} 进行累加
        for prot_id, sp_vector in protein_vectors.items():
            for subtype in subtypes:
                if subtype in sp_vector:
                    s_locus[subtype] += sp_vector[subtype]
                    
        # 寻找最高分 (argmax S_Locus)
        best_subtype = max(s_locus, key=s_locus.get)
        best_score = s_locus[best_subtype]
        
        if best_score > 0:
            # 简化版得分转换，可根据实际需求调整阈值映射
            adaptation_score = min(1.0, best_score / 100.0) 
            interference_score = adaptation_score
            return {'adaptation_score': adaptation_score, 'interference_score': interference_score}, best_subtype
            
        return {'adaptation_score': 0.0, 'interference_score': 0.0}, ""
    
    def _load_rules(self, rules_path: str) -> Dict:
        """加载规则文件"""
        try:
            with open(rules_path, 'r') as f:
                return json.load(f)
        except Exception as e:
            print(f"Failed to load rules from {rules_path}: {e}")
            # 返回默认规则
            return self._get_default_rules()
    
    def _load_cas_scoring_matrix(self, scoring_path: str) -> pd.DataFrame:
        """加载专家规则权重矩阵"""
        try:
            if not Path(scoring_path).exists():
                print(f"CasScoring matrix not found at {scoring_path}, using simplified rules")
                return None
            
            # 读取CSV文件，第一列为索引（Hmm列）
            df = pd.read_csv(scoring_path, index_col=0, na_values=['', ' ', 'NaN'])
            # 将空值替换为0
            df = df.fillna(0.0)
            print(f"✅ Loaded CasScoring matrix with shape: {df.shape}")
            return df
        except Exception as e:
            print(f"Failed to load CasScoring matrix from {scoring_path}: {e}")
            return None
    
    def _get_default_rules(self) -> Dict:
        """返回默认规则"""
        return {
            "cas_families": ["Cas1", "Cas2", "Cas3", "Cas4", "Cas5", "Cas6", "Cas7", "Cas8", "Cas9", "Cas10", "Cas12", "Cas13"],
            "subtypes": [f"Type_{i}" for i in range(1, 34)],
            "rules": {}  # 简化的规则矩阵
        }
    
    def classify_system(self, 
                       repeat_subtype: str, 
                       cas_families: List[str],
                       hmm_results: List[Dict],
                       protein_vectors: Dict[str, Dict[str, float]] = None) -> Dict:
        """基于规则对CRISPR-Cas系统进行分类"""
        
        if self.cas_scoring_matrix is not None and protein_vectors:
            expert_score, best_subtype = self._aggregate_locus_score(protein_vectors)
            adaptation_score = expert_score.get('adaptation_score', 0.0)
            interference_score = expert_score.get('interference_score', 0.0)
            final_subtype = best_subtype if best_subtype else repeat_subtype
        else:
            adaptation_score = self._calculate_adaptation_score(cas_families)
            interference_score = self._calculate_interference_score(cas_families, hmm_results)
            final_subtype = repeat_subtype
        
        # 结合重复序列预测
        repeat_confidence = 0.8  # 假设重复序列分类的置信度
        
        # 综合评分
        final_score = (adaptation_score + interference_score + repeat_confidence) / 3.0
        
        # 确定最终亚型
        if final_score >= self.confidence_threshold:
            predicted_subtype = final_subtype
        else:
            predicted_subtype = "Unclassified"
        
        return {
            'predicted_subtype': predicted_subtype,
            'confidence': final_score,
            'adaptation_score': adaptation_score,
            'interference_score': interference_score,
            'repeat_subtype': repeat_subtype,
            'detected_cas_families': cas_families,
            'hmm_matches': [r['hmm_family'] for r in hmm_results]
        }
    
    def _calculate_expert_score(self, cas_families: List[str], hmm_results: List[Dict]) -> Tuple[Dict, str]:
        """使用专家规则权重矩阵计算评分"""
        if self.cas_scoring_matrix is None:
            return {'adaptation_score': 0.0, 'interference_score': 0.0}, ""
        
        # 获取所有可能的亚型列表
        subtypes = list(self.cas_scoring_matrix.columns)
        
        # 合并cas_families和hmm_results中的家族信息
        all_detected_families = set(cas_families)
        for result in hmm_results:
            if 'hmm_family' in result and result['hmm_family']:
                all_detected_families.add(result['hmm_family'])
        
        # 计算每个亚型的总分
        subtype_scores = {}
        for subtype in subtypes:
            total_score = 0.0
            for family in all_detected_families:
                # 在权重矩阵中查找该家族对该亚型的评分
                if family in self.cas_scoring_matrix.index:
                    score = self.cas_scoring_matrix.loc[family, subtype]
                    if pd.notna(score) and score != 0:
                        total_score += float(score)
            subtype_scores[subtype] = total_score
        
        # 找到最高分的亚型
        if subtype_scores:
            best_subtype = max(subtype_scores, key=subtype_scores.get)
            best_score = subtype_scores[best_subtype]
            
            # 根据专家矩阵计算适应和干扰模块得分
            # 这里简化处理：将总分按比例分配到两个模块
            # 实际应用中可以根据具体的亚型类型来分配
            adaptation_score = min(1.0, max(0.0, best_score / 10.0))  # 归一化到0-1范围
            interference_score = adaptation_score
            
            return {'adaptation_score': adaptation_score, 'interference_score': interference_score}, best_subtype
        
        return {'adaptation_score': 0.0, 'interference_score': 0.0}, ""
    
    def _calculate_adaptation_score(self, cas_families: List[str]) -> float:
        """计算适应模块得分"""
        # 适应模块通常需要Cas1和Cas2
        required_adaptation = {"Cas1", "Cas2"}
        detected_set = set(cas_families)
        
        if required_adaptation.issubset(detected_set):
            return 1.0
        elif "Cas1" in detected_set:
            return 0.7
        elif "Cas2" in detected_set:
            return 0.5
        else:
            return 0.2
    
    def _calculate_interference_score(self, cas_families: List[str], hmm_results: List[Dict]) -> float:
        """计算干扰模块得分"""
        # 干扰模块需要特定的效应蛋白
        interference_proteins = {"Cas3", "Cas9", "Cas10", "Cas12", "Cas13"}
        detected_set = set(cas_families)
        
        # 检查是否有干扰模块蛋白
        interference_detected = detected_set.intersection(interference_proteins)
        
        if len(interference_detected) > 0:
            return 1.0
        elif len(hmm_results) > 0:
            # 如果HMM匹配到相关家族，给予部分分数
            return 0.6
        else:
            return 0.1