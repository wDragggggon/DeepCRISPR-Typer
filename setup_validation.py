#!/usr/bin/env python3
"""
Setup validation script for DeepCRISPR-Typer.
This script helps users verify their installation and identify missing components.
"""

import os
import sys
import shutil
from pathlib import Path

def check_python_dependencies():
    """Check if required Python packages are installed"""
    required_packages = [
        'torch', 'biopython', 'pandas', 'numpy', 'sklearn', 
        'matplotlib', 'seaborn', 'xgboost', 'fair_esm'
    ]
    
    print("🔍 Checking Python dependencies...")
    missing_packages = []
    
    for package in required_packages:
        try:
            __import__(package)
            print(f"  ✅ {package}")
        except ImportError:
            missing_packages.append(package)
            print(f"  ❌ {package} (missing)")
    
    return missing_packages

def check_external_tools():
    """Check if external tools are available"""
    from config.paths_config import MINCED_PATH, PRODIGAL_PATH, HMMSCAN_PATH
    
    print("\n🔍 Checking external tools...")
    missing_tools = []
    
    # Check MinCED
    if os.path.exists(MINCED_PATH):
        print(f"  ✅ MinCED: {MINCED_PATH}")
    else:
        missing_tools.append(f"MinCED ({MINCED_PATH})")
        print(f"  ❌ MinCED: {MINCED_PATH} (not found)")
    
    # Check Prodigal
    if PRODIGAL_PATH == "prodigal":
        # Try to find in PATH
        prodigal_path = shutil.which("prodigal")
        if prodigal_path:
            print(f"  ✅ Prodigal: {prodigal_path}")
        else:
            missing_tools.append("Prodigal (not in PATH)")
            print("  ❌ Prodigal: not found in system PATH")
    else:
        if os.path.exists(PRODIGAL_PATH):
            print(f"  ✅ Prodigal: {PRODIGAL_PATH}")
        else:
            missing_tools.append(f"Prodigal ({PRODIGAL_PATH})")
            print(f"  ❌ Prodigal: {PRODIGAL_PATH} (not found)")
    
    # Check HMMER
    if HMMSCAN_PATH == "hmmscan":
        hmmscan_path = shutil.which("hmmscan")
        if hmmscan_path:
            print(f"  ✅ HMMER: {hmmscan_path}")
        else:
            missing_tools.append("HMMER (not in PATH)")
            print("  ❌ HMMER: not found in system PATH")
    else:
        if os.path.exists(HMMSCAN_PATH):
            print(f"  ✅ HMMER: {HMMSCAN_PATH}")
        else:
            missing_tools.append(f"HMMER ({HMMSCAN_PATH})")
            print(f"  ❌ HMMER: {HMMSCAN_PATH} (not found)")
    
    return missing_tools

def check_data_files():
    """Check if required data files and models exist"""
    from config.paths_config import (
        CNN_ATT_MODEL_PATH, TEMC_CAS_MODEL_PATH, HMM_DB_PATH,
        ADAPTATION_RULES_PATH, INTERFERENCE_RULES_PATH, CAS_SCORING_PATH
    )
    
    print("\n🔍 Checking data files and models...")
    missing_files = []
    
    # Required files
    required_files = {
        "CNN-Att Model": CNN_ATT_MODEL_PATH,
        "HMM Database": HMM_DB_PATH,
        "Adaptation Rules": ADAPTATION_RULES_PATH,
        "Interference Rules": INTERFERENCE_RULES_PATH,
        "CasScoring Matrix": CAS_SCORING_PATH
    }
    
    for name, path in required_files.items():
        if os.path.exists(path):
            print(f"  ✅ {name}: {path}")
        else:
            missing_files.append(f"{name} ({path})")
            print(f"  ❌ {name}: {path} (not found)")
    
    # Optional files
    optional_files = {
        "TEMC-Cas Model": TEMC_CAS_MODEL_PATH
    }
    
    for name, path in optional_files.items():
        if os.path.exists(path):
            print(f"  ✅ {name}: {path}")
        else:
            print(f"  ⚠️  {name}: {path} (optional, not found)")
    
    return missing_files

def main():
    """Main validation function"""
    print("🚀 DeepCRISPR-Typer Setup Validation")
    print("=" * 50)
    
    # Add project to Python path
    project_root = Path(__file__).parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    
    try:
        missing_deps = check_python_dependencies()
        missing_tools = check_external_tools()
        missing_data = check_data_files()
        
        print("\n" + "=" * 50)
        print("📋 SUMMARY")
        
        if missing_deps:
            print(f"\n❌ Missing Python packages ({len(missing_deps)}):")
            for pkg in missing_deps:
                print(f"   - {pkg}")
        
        if missing_tools:
            print(f"\n❌ Missing external tools ({len(missing_tools)}):")
            for tool in missing_tools:
                print(f"   - {tool}")
        
        if missing_data:
            print(f"\n❌ Missing data files/models ({len(missing_data)}):")
            for file in missing_data:
                print(f"   - {file}")
        
        if not (missing_deps or missing_tools or missing_data):
            print("\n✅ All required components are properly installed!")
            print("🎉 You're ready to use DeepCRISPR-Typer!")
        else:
            print(f"\n⚠️  Please install the missing components before proceeding.")
            print("📚 Refer to SETUP.md for detailed installation instructions.")
    
    except Exception as e:
        print(f"\n💥 Error during validation: {e}")
        print("🔍 Please check your Python environment and configuration.")

if __name__ == "__main__":
    main()