from __future__ import annotations
from dataclasses import dataclass
import hashlib
from typing import Dict, List, Optional

from factorforge.analysis.metrics import load_codon_usage_table
from factorforge.constraints.type_iis import get_canonical_forbidden_motifs
from factorforge.discovery.filter import HardConstraintConfig, HardConstraintFilter
from factorforge.engines.sllm.automaton import AutomatonCompiler
from factorforge.engines.sllm.interfaces import (
    AutomatonConstraintProcessor,
    ConstraintProcessorPipeline,
    HomopolymerConstraintProcessor,
    SynonymousMaskProcessor,
)


@dataclass(frozen=True)
class ResolvedDesignContract:
    """Immutable, mathematically unified DesignContract snapshot shared by sLLM and DP engines."""

    contract_id: str
    contract_digest: str
    host: str
    hard_config: HardConstraintConfig
    decode_pipeline: ConstraintProcessorPipeline
    final_validator: HardConstraintFilter
    codon_weights: Dict[str, float]
    expected_stop: str = "TAA"
    upstream_flank: str = ""
    downstream_flank: str = ""


class DesignContractCompiler:
    """Compiles canonical HardConstraintConfig into executable sLLM constraint processor pipelines and resolved contracts."""

    @classmethod
    def compile_from_config(
        cls,
        config: Optional[HardConstraintConfig] = None,
        custom_motifs: Optional[List[str]] = None,
    ) -> ConstraintProcessorPipeline:
        """Builds a deterministic ConstraintProcessorPipeline from standard constraint configuration."""
        cfg = config or HardConstraintConfig()

        # 1. Exact Amino Acid translation processor (Hard Veto)
        syn_proc = SynonymousMaskProcessor()

        # 2. Compile canonical forbidden motifs (Type IIS + custom) via Aho-Corasick (Hard Veto)
        canonical_motifs = get_canonical_forbidden_motifs(cfg.forbidden_type_iis_enzymes)
        all_motifs = set(canonical_motifs)

        if cfg.forbidden_other_motifs:
            for m in cfg.forbidden_other_motifs:
                all_motifs.add(m.strip().upper())

        if custom_motifs:
            for m in custom_motifs:
                all_motifs.add(m.strip().upper())

        # Compile Aho-Corasick DFA automaton with reverse-complement support
        automaton = AutomatonCompiler.compile(all_motifs, include_rc=True)
        auto_proc = AutomatonConstraintProcessor(automaton)

        # 3. Homopolymer limit processor (Hard Veto)
        homo_proc = HomopolymerConstraintProcessor(max_run=cfg.homopolymer_max_run)

        return ConstraintProcessorPipeline(
            [
                syn_proc,
                auto_proc,
                homo_proc,
            ]
        )

    @classmethod
    def resolve_contract(
        cls,
        host: str = "nbenthamiana",
        config: Optional[HardConstraintConfig] = None,
        custom_motifs: Optional[List[str]] = None,
        expected_stop: str = "TAA",
        upstream_flank: str = "",
        downstream_flank: str = "",
        contract_id: Optional[str] = None,
    ) -> ResolvedDesignContract:
        """Constructs an immutable ResolvedDesignContract binding decode pipeline, DP parameters, and final validator."""
        cfg = config or HardConstraintConfig()
        pipeline = cls.compile_from_config(cfg, custom_motifs=custom_motifs)
        validator = HardConstraintFilter(cfg)

        # Load codon weights for host
        table = load_codon_usage_table()
        codon_weights = dict(table.codon_weights)

        cid = contract_id or f"DC-{host}-{expected_stop}"
        digest_payload = (
            f"{host}:{cfg.global_gc_min}:{cfg.global_gc_max}:{cfg.homopolymer_max_run}:"
            f"{sorted(list(cfg.forbidden_type_iis_enzymes))}:{sorted(cfg.forbidden_other_motifs)}:"
            f"{expected_stop}:{upstream_flank}:{downstream_flank}"
        )
        digest = f"sha256:{hashlib.sha256(digest_payload.encode('utf-8')).hexdigest()}"

        return ResolvedDesignContract(
            contract_id=cid,
            contract_digest=digest,
            host=host,
            hard_config=cfg,
            decode_pipeline=pipeline,
            final_validator=validator,
            codon_weights=codon_weights,
            expected_stop=expected_stop,
            upstream_flank=upstream_flank,
            downstream_flank=downstream_flank,
        )
