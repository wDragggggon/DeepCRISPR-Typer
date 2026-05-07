import argparse
import os
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent / "src"))

from src.pipeline import DeepCRISPRTyperPipeline
from config.paths_config import validate_paths


def main():
    parser = argparse.ArgumentParser(
        description="DeepCRISPR-Typer: CRISPR-Cas System Classification Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python main.py -i input_genome.fasta
  python main.py -i input_genome.fasta -o my_results
        """
    )

    parser.add_argument(
        "-i", "--input",
        help="Input genome FASTA file path"
    )

    parser.add_argument(
        "-o", "--output",
        help="Output prefix (default: use input filename)"
    )

    parser.add_argument(
        "--validate",
        action="store_true",
        help="Validate environment configuration only without running the pipeline"
    )

    args = parser.parse_args()

    # Validation mode does not require input file
    if args.validate:
        if validate_paths():
            print("✅ Environment configuration validation passed!")
        else:
            print("❌ Environment configuration validation failed. Please check required paths and tools.")
            sys.exit(1)
        return

    # Normal mode requires input file
    if not args.input:
        parser.error("the following arguments are required: -i/--input")

    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input FASTA file not found: {args.input}")

    # Run pipeline
    try:
        pipeline = DeepCRISPRTyperPipeline()
        results = pipeline.run(args.input, args.output)

        if results:
            print(f"\n🎉 Successfully processed {len(results)} CRISPR arrays!")
        else:
            print("\n⚠️ No CRISPR arrays detected or processing failed.")

    except Exception as e:
        print(f"❌ Pipeline execution failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()