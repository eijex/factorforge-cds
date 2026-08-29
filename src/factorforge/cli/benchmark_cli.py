import argparse
import json
import os
from factorforge.benchmark.models import BenchmarkRunConfig, BenchmarkTarget
from factorforge.benchmark.runner import BenchmarkRunner

def generate_markdown_report(report_json: dict, output_path: str):
    """Generate a human-readable Markdown report from the JSON artifact."""
    lines = [
        f"# Benchmark Suite Report: {report_json['suite_name']}",
        f"**Run Timestamp**: {report_json['run_timestamp']}",
        "",
        "## Summary (Pass Rate)",
        "*(Format: Passed / Total Targets)*"
    ]
    
    for engine, rate in report_json['pass_counts'].items():
        lines.append(f"- **{engine}**: {rate}")
        
    lines.append("")
    lines.append("## Average Metrics")
    lines.append("| Engine | Average CAI | Average GC% |")
    lines.append("|--------|-------------|-------------|")
    
    for engine in report_json['average_cai']:
        cai = report_json['average_cai'][engine]
        gc = report_json['average_gc'][engine] * 100
        lines.append(f"| {engine} | {cai:.3f} | {gc:.1f}% |")
        
    lines.append("")
    lines.append("## Target Results")
    
    for res in report_json['results']:
        engine = res['engine_name']
        target = res['target_id']
        passed = not any(
            c['enforcement'] in ('hard_fail', 'gate') and c['result'] == 'fail'
            for c in res['evaluation']['checks']
        )
        lines.append(f"### Target: {target} | Engine: {engine} | Passed: {'✅' if passed else '❌'}")
        lines.append(f"- **Runtime**: {res['runtime_seconds']:.2f}s")
        lines.append(f"- **GC%**: {res['evaluation']['metrics']['gc_percent'] * 100:.1f}%")
        if not passed:
            lines.append("- **Failing Checks**:")
            for check in res['evaluation']['checks']:
                if check['result'] != "PASS":
                    lines.append(f"  - {check['check_name']} ({check['enforcement']}): {check['message']}")
        lines.append("")

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description="FactorForge Benchmark Harness")
    parser.add_argument("--engines", type=str, required=True, help="Comma-separated list of engines (e.g. profile,dp,lm)")
    parser.add_argument("--suite", type=str, default="Engineering Benchmark Suite")
    parser.add_argument("--host", type=str, default="nbenthamiana")
    parser.add_argument("--outdir", type=str, default="reports/benchmarks")
    args = parser.parse_args()

    engine_names = [e.strip() for e in args.engines.split(",")]
    
    # Mock Reference Targets for Engineering Suite
    targets = [
        BenchmarkTarget(target_id="eGFP", protein_sequence="MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"),
        BenchmarkTarget(target_id="Humira_HC", protein_sequence="EVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK")
    ]

    config = BenchmarkRunConfig(
        suite_name=args.suite,
        host=args.host,
        target_gc_min=0.40,
        target_gc_max=0.47,
        terminal_stop_policy="append",
        forbidden_type_iis=["BsaI", "BsmBI", "BpiI"]
    )

    runner = BenchmarkRunner(config)
    report = runner.run_suite(targets, engine_names)
    
    os.makedirs(args.outdir, exist_ok=True)
    report_json = report.model_dump()
    
    json_path = os.path.join(args.outdir, "benchmark_report.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(report_json, f, indent=2)
        
    md_path = os.path.join(args.outdir, "benchmark_report.md")
    generate_markdown_report(report_json, md_path)
    
    print(f"Benchmark complete. Machine-readable artifact saved to {json_path}. Human-readable report saved to {md_path}.")

if __name__ == "__main__":
    main()
