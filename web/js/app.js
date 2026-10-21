/**
 * FactorForge Application Logic
 * Vanilla JS for the FactorForge CDS design interface.
 */

const API_ENDPOINT = '/api/optimize';
const ENABLE_MOCK = window.FACTORFORGE_ENABLE_MOCK === true;
// Maps internal engine table names (incl. HOST_MAP aliases like 'ntabacum')
// back to a human label for results and legacy local-history display. The web
// design workflow itself is intentionally fixed to N. benthamiana.
const HOST_LABELS = {
    nbenthamiana: 'N. benthamiana',
    by2: 'Tobacco BY-2',
    ntabacum: 'Tobacco BY-2',
    human_stemcell: 'Human Stem Cell (hMSC / iPSC)'
};
// GC reference band per host. Populated live from GET /api/optimize's
// host_metadata[host].gc_range (api/optimize.py's _default_gc_constraints /
// resolve_host_gc_range — the single source of truth, registry-synced).
// OFFLINE_GC_RANGES is only a same-values fallback for offline/dev mode when
// the API is unreachable (v3.3.0) — never the other
// way around, so this file must not be the place a future band change is made.
const OFFLINE_GC_RANGES = {
    nbenthamiana: { gc_min: 40.0, gc_max: 47.0 },
    by2: { gc_min: 55.0, gc_max: 65.0 },
    ntabacum: { gc_min: 55.0, gc_max: 65.0 },
    human_stemcell: { gc_min: 45.0, gc_max: 50.0 }
};
// Recognition sequences mirror Domesticator.ASSEMBLY_STANDARDS["golden_gate"]
// and design_review.TYPE_IIS_SITES["SapI"]. scan_rc lets the API check both strands.
const TYPE_IIS_PRESETS = Object.freeze({
    BsaI: 'GGTCTC',
    BpiI: 'GAAGAC',
    BsmBI: 'CGTCTC',
    SapI: 'GAAGAGC'
});
let hostGcRanges = {};
let hostMetadata = {};
let apiCapabilities = {};

const CODON_TABLE = Object.freeze({
    TTT:'F',TTC:'F',TTA:'L',TTG:'L',TCT:'S',TCC:'S',TCA:'S',TCG:'S',TAT:'Y',TAC:'Y',TAA:'*',TAG:'*',TGT:'C',TGC:'C',TGA:'*',TGG:'W',
    CTT:'L',CTC:'L',CTA:'L',CTG:'L',CCT:'P',CCC:'P',CCA:'P',CCG:'P',CAT:'H',CAC:'H',CAA:'Q',CAG:'Q',CGT:'R',CGC:'R',CGA:'R',CGG:'R',
    ATT:'I',ATC:'I',ATA:'I',ATG:'M',ACT:'T',ACC:'T',ACA:'T',ACG:'T',AAT:'N',AAC:'N',AAA:'K',AAG:'K',AGT:'S',AGC:'S',AGA:'R',AGG:'R',
    GTT:'V',GTC:'V',GTA:'V',GTG:'V',GCT:'A',GCC:'A',GCA:'A',GCG:'A',GAT:'D',GAC:'D',GAA:'E',GAG:'E',GGT:'G',GGC:'G',GGA:'G',GGG:'G'
});
const AA_NAMES = Object.freeze({A:'Ala',R:'Arg',N:'Asn',D:'Asp',C:'Cys',Q:'Gln',E:'Glu',G:'Gly',H:'His',I:'Ile',L:'Leu',K:'Lys',M:'Met',F:'Phe',P:'Pro',S:'Ser',T:'Thr',W:'Trp',Y:'Tyr',V:'Val','*':'Stop'});

function getGcRange(hostId) {
    return hostGcRanges[hostId] || OFFLINE_GC_RANGES[hostId] || OFFLINE_GC_RANGES.nbenthamiana;
}

let validationRegistry = [];
const HISTORY_SCHEMA_VERSION = 3;
// Custom SOP management is intentionally paused. Retained code can be re-enabled later.
const ENABLE_CUSTOM_SOP = false;
const SOP_STORAGE_KEY = 'factorforge_active_sop_v2';
const SOP_SCHEMA = 'factorforge-sop-v1';
const SOP_MODES = Object.freeze(['required', 'preferred', 'ignored']);
const SOP_PRESETS = Object.freeze({
    preset_a: {
        id: 'preset_a',
        name: 'Preset A · Conservative Baseline (Default)',
        tag: 'Preset A Active',
        file: '/examples/factorforge-conservative-sop.yaml',
        description: 'Pre-synthesis safety baseline. Enforces 5′ initiation ramp (nt 1–45 low-GC clamp), homopolymer veto (≤5 nt), and standard Type IIS avoidance.',
        profile: {
            $schema: SOP_SCHEMA,
            profile_id: 'preset_a_conservative_baseline',
            sop_name: 'Preset A · Conservative Baseline (Default)',
            version: '1.2.0',
            author: 'Eijex',
            status: 'STANDARD_TEMPLATE',
            derived_from: 'FactorForge conservative pre-synthesis baseline; precautionary plant expression heuristics',
            scope_note: 'Conservative pre-synthesis baseline. Custom laboratory SOP upload or configuration is recommended.',
            default_enforcement: 'IGNORE',
            unknown_rule_policy: 'ERROR',
            rules: {
                'assembly.type_iis.bsai.v1': 'HARD_FAIL',
                'assembly.type_iis.bpii.v1': 'HARD_FAIL',
                'assembly.type_iis.bsmbi.v1': 'HARD_FAIL',
                'biological.reading_frame.v1': 'HARD_FAIL',
                'rna.cryptic_splice.v1': 'WARNING',
                'rna.polya_motifs.v1': 'WARNING',
                'rna.au_rich_elements.v1': 'WARNING',
                'synthesis.gc_extremes.v1': 'WARNING',
                'policy.synthesis.homopolymer.v1': 'WARNING',
                'policy.initiation_mfe.v1': 'WARNING'
            },
            workflow: {
                design_method: 'dp_v2_1_1',
                codon_policy: 'plantform_balanced_codon_v1',
                watermark: false,
                comparison_methods: ['high_cai', 'gc_target', 'assembly_friendly'],
                sequence_requirements: {
                    reproducibility_seed: null,
                    kozak_optimization: false,
                    suppress_tpa: false,
                    type_iis_enzymes: ['BsaI', 'BpiI', 'BsmBI'],
                    custom_restriction_sites: [],
                    moclo_template: false
                },
                review_policy: {
                    cai: 'PREFERRED', overall_gc: 'PREFERRED', local_gc: 'PREFERRED',
                    type_iis: 'REQUIRED', repeats: 'PREFERRED', homopolymers: 'PREFERRED',
                    forbidden_motifs: 'PREFERRED'
                }
            }
        }
    },
    preset_b: {
        id: 'preset_b',
        name: 'Preset B · Therapeutic Biologics (Antibodies & Complex Proteins)',
        tag: 'Preset B Active',
        file: '/examples/factorforge-sop-preset-b-therapeutic-biologics.yaml',
        description: 'Biologics expression profile. Enforces balanced non-CpG codons (Arg 50/50, Ser 5-way), plant Kozak initiation, mRNA stability (TpA suppression), and comprehensive Type IIS avoidance.',
        profile: {
            $schema: SOP_SCHEMA,
            profile_id: 'preset_b_therapeutic_biologics',
            sop_name: 'Preset B · Therapeutic Biologics (Antibodies & Complex Proteins)',
            version: '1.2.0',
            author: 'Eijex',
            status: 'PRESET_TEMPLATE',
            derived_from: 'Plant expression heuristics for therapeutic antibodies and complex recombinant glycoproteins',
            scope_note: 'Optimized sequence policy for therapeutic antibodies and complex proteins in N. benthamiana',
            default_enforcement: 'IGNORE',
            unknown_rule_policy: 'ERROR',
            rules: {
                'assembly.type_iis.bsai.v1': 'HARD_FAIL',
                'assembly.type_iis.bpii.v1': 'HARD_FAIL',
                'assembly.type_iis.bsmbi.v1': 'HARD_FAIL',
                'assembly.type_iis.sapi.v1': 'HARD_FAIL',
                'biological.reading_frame.v1': 'HARD_FAIL',
                'rna.cryptic_splice.v1': 'HARD_FAIL',
                'rna.polya_motifs.v1': 'HARD_FAIL',
                'rna.au_rich_elements.v1': 'WARNING',
                'synthesis.gc_extremes.v1': 'WARNING',
                'policy.synthesis.homopolymer.v1': 'WARNING',
                'policy.initiation_mfe.v1': 'WARNING'
            },
            workflow: {
                design_method: 'dp_v2_1_1',
                codon_policy: 'plantform_balanced_codon_v1',
                watermark: false,
                comparison_methods: ['high_cai', 'gc_target', 'assembly_friendly'],
                sequence_requirements: {
                    reproducibility_seed: null,
                    kozak_optimization: true,
                    suppress_tpa: true,
                    type_iis_enzymes: ['BsaI', 'BpiI', 'BsmBI', 'SapI'],
                    custom_restriction_sites: [],
                    moclo_template: false
                },
                review_policy: {
                    cai: 'REQUIRED', overall_gc: 'REQUIRED', local_gc: 'REQUIRED',
                    type_iis: 'REQUIRED', repeats: 'PREFERRED', homopolymers: 'PREFERRED',
                    forbidden_motifs: 'REQUIRED'
                }
            }
        }
    },
    preset_c: {
        id: 'preset_c',
        name: 'Preset C · MoClo & GoldenBraid Modular Assembly',
        tag: 'Preset C Active',
        file: '/examples/factorforge-sop-preset-c-moclo-goldenbraid.yaml',
        description: 'MoClo Level 0 & GoldenBraid 2.0 standards. Complete avoidance of BsaI/BpiI/BsmBI/SapI, MoClo overhang compatibility (5′ AATG, 3′ GCTT), and strict repeat filtering.',
        profile: {
            $schema: SOP_SCHEMA,
            profile_id: 'preset_c_moclo_goldenbraid',
            sop_name: 'Preset C · MoClo & GoldenBraid Modular Assembly',
            version: '1.2.0',
            author: 'Eijex',
            status: 'PRESET_TEMPLATE',
            derived_from: 'Modular Cloning (MoClo) and GoldenBraid 2.0 plant assembly standards',
            scope_note: 'Standard sequence policy for multi-gene and modular plant vector cloning',
            default_enforcement: 'IGNORE',
            unknown_rule_policy: 'ERROR',
            rules: {
                'assembly.type_iis.bsai.v1': 'HARD_FAIL',
                'assembly.type_iis.bpii.v1': 'HARD_FAIL',
                'assembly.type_iis.bsmbi.v1': 'HARD_FAIL',
                'assembly.type_iis.sapi.v1': 'HARD_FAIL',
                'biological.reading_frame.v1': 'HARD_FAIL',
                'rna.cryptic_splice.v1': 'WARNING',
                'rna.polya_motifs.v1': 'WARNING',
                'rna.au_rich_elements.v1': 'WARNING',
                'synthesis.gc_extremes.v1': 'WARNING',
                'policy.synthesis.homopolymer.v1': 'WARNING',
                'policy.initiation_mfe.v1': 'WARNING'
            },
            workflow: {
                design_method: 'assembly_friendly',
                codon_policy: 'plantform_balanced_codon_v1',
                watermark: false,
                comparison_methods: ['dp_v2_1_1', 'high_cai', 'gc_target'],
                sequence_requirements: {
                    reproducibility_seed: null,
                    kozak_optimization: false,
                    suppress_tpa: false,
                    type_iis_enzymes: ['BsaI', 'BpiI', 'BsmBI', 'SapI'],
                    custom_restriction_sites: [],
                    moclo_template: true
                },
                review_policy: {
                    cai: 'PREFERRED', overall_gc: 'PREFERRED', local_gc: 'PREFERRED',
                    type_iis: 'REQUIRED', repeats: 'REQUIRED', homopolymers: 'REQUIRED',
                    forbidden_motifs: 'REQUIRED'
                }
            }
        }
    },
    preset_d: {
        id: 'preset_d',
        name: 'Preset D · High-Yield Agroinfiltration (Transient Overexpression)',
        tag: 'Preset D Active',
        file: '/examples/factorforge-sop-preset-d-high-yield-agro.yaml',
        description: 'High-yield transient expression via agroinfiltration. Max CAI codon adaptation, optimized Kozak initiation, and mRNA destabilization suppression.',
        profile: {
            $schema: SOP_SCHEMA,
            profile_id: 'preset_d_high_yield_agro',
            sop_name: 'Preset D · High-Yield Agroinfiltration Overexpression',
            version: '1.2.0',
            author: 'Eijex',
            status: 'PRESET_TEMPLATE',
            derived_from: 'Plant transient overexpression high-throughput translation heuristics',
            scope_note: 'Maximum yield protein expression policy for N. benthamiana leaf agroinfiltration',
            default_enforcement: 'IGNORE',
            unknown_rule_policy: 'ERROR',
            rules: {
                'assembly.type_iis.bsai.v1': 'HARD_FAIL',
                'assembly.type_iis.bpii.v1': 'HARD_FAIL',
                'biological.reading_frame.v1': 'HARD_FAIL',
                'rna.cryptic_splice.v1': 'WARNING',
                'rna.polya_motifs.v1': 'WARNING',
                'rna.au_rich_elements.v1': 'WARNING',
                'synthesis.gc_extremes.v1': 'WARNING',
                'policy.synthesis.homopolymer.v1': 'WARNING',
                'policy.initiation_mfe.v1': 'WARNING'
            },
            workflow: {
                design_method: 'high_cai',
                codon_policy: 'max_cai',
                watermark: false,
                comparison_methods: ['dp_v2_1_1', 'gc_target', 'assembly_friendly'],
                sequence_requirements: {
                    reproducibility_seed: null,
                    kozak_optimization: true,
                    suppress_tpa: true,
                    type_iis_enzymes: ['BsaI', 'BpiI'],
                    custom_restriction_sites: [],
                    moclo_template: false
                },
                review_policy: {
                    cai: 'REQUIRED', overall_gc: 'PREFERRED', local_gc: 'PREFERRED',
                    type_iis: 'REQUIRED', repeats: 'PREFERRED', homopolymers: 'PREFERRED',
                    forbidden_motifs: 'PREFERRED'
                }
            }
        }
    },
    preset_e: {
        id: 'preset_e',
        name: 'Preset E · Synthesis House Feasibility (Twist / IDT Pass-Through)',
        tag: 'Preset E Active',
        file: '/examples/factorforge-sop-preset-e-synthesis-feasibility.yaml',
        description: 'Optimized for commercial DNA synthesis (Twist, IDT, GenScript) complexity metrics. Distributes host codon frequencies and strictly limits GC extremes and homopolymers.',
        profile: {
            $schema: SOP_SCHEMA,
            profile_id: 'preset_e_synthesis_feasibility',
            sop_name: 'Preset E · Synthesis House Feasibility (Twist / IDT)',
            version: '1.2.0',
            author: 'Eijex',
            status: 'PRESET_TEMPLATE',
            derived_from: 'Commercial DNA synthesis manufacturing complexity guidelines (Twist / IDT / GenScript)',
            scope_note: 'Sequence policy targeting commercial synthesis approval without complexity delays',
            default_enforcement: 'IGNORE',
            unknown_rule_policy: 'ERROR',
            rules: {
                'assembly.type_iis.bsai.v1': 'HARD_FAIL',
                'assembly.type_iis.bsmbi.v1': 'HARD_FAIL',
                'biological.reading_frame.v1': 'HARD_FAIL',
                'rna.cryptic_splice.v1': 'WARNING',
                'rna.polya_motifs.v1': 'WARNING',
                'rna.au_rich_elements.v1': 'WARNING',
                'synthesis.gc_extremes.v1': 'HARD_FAIL',
                'policy.synthesis.homopolymer.v1': 'HARD_FAIL',
                'policy.initiation_mfe.v1': 'WARNING'
            },
            workflow: {
                design_method: 'feasibility_best',
                codon_policy: 'host_frequency',
                watermark: false,
                comparison_methods: ['dp_v2_1_1', 'high_cai', 'assembly_friendly'],
                sequence_requirements: {
                    reproducibility_seed: null,
                    kozak_optimization: false,
                    suppress_tpa: false,
                    type_iis_enzymes: ['BsaI', 'BsmBI'],
                    custom_restriction_sites: [],
                    moclo_template: false
                },
                review_policy: {
                    cai: 'PREFERRED', overall_gc: 'REQUIRED', local_gc: 'REQUIRED',
                    type_iis: 'REQUIRED', repeats: 'REQUIRED', homopolymers: 'REQUIRED',
                    forbidden_motifs: 'PREFERRED'
                }
            }
        }
    }
});

const DEFAULT_SOP_PROFILE = Object.freeze(SOP_PRESETS.preset_a.profile);

function cloneJson(value) {
    return JSON.parse(JSON.stringify(value));
}

function parseYamlScalar(raw) {
    const value = raw.trim();
    if (!value || value === 'null' || value === '~') return null;
    if (/^(true|yes)$/i.test(value)) return true;
    if (/^(false|no)$/i.test(value)) return false;
    if (/^-?(?:\d+\.?\d*|\.\d+)$/.test(value)) return Number(value);
    if (value.startsWith('[') && value.endsWith(']')) {
        const inner = value.slice(1, -1).trim();
        return inner ? inner.split(',').map(item => parseYamlScalar(item)) : [];
    }
    if ((value.startsWith('"') && value.endsWith('"')) || (value.startsWith("'") && value.endsWith("'"))) {
        return value.startsWith('"') ? JSON.parse(value) : value.slice(1, -1).replace(/''/g, "'");
    }
    return value;
}

function stripYamlComment(line) {
    let quote = null;
    for (let index = 0; index < line.length; index += 1) {
        const char = line[index];
        if ((char === '"' || char === "'") && line[index - 1] !== '\\') quote = quote === char ? null : (quote || char);
        if (char === '#' && !quote && (index === 0 || /\s/.test(line[index - 1]))) return line.slice(0, index);
    }
    return line;
}

function parseSopYaml(text) {
    if (/\t/.test(text)) throw new Error('YAML indentation must use spaces, not tabs.');
    if (/(^|\s)[&*!][^\s]*/m.test(text) || /^\s*[>|]/m.test(text)) throw new Error('YAML anchors, tags, aliases, and multiline blocks are not supported.');
    const root = {};
    const stack = [{ indent: -1, value: root }];
    text.split(/\r?\n/).forEach((source, lineIndex) => {
        const uncommented = stripYamlComment(source).replace(/\s+$/, '');
        if (!uncommented.trim() || uncommented.trim() === '---') return;
        const indent = uncommented.length - uncommented.trimStart().length;
        const content = uncommented.trimStart();
        const separator = content.indexOf(':');
        if (separator < 1) throw new Error(`Invalid YAML on line ${lineIndex + 1}.`);
        const key = content.slice(0, separator).trim().replace(/^['"]|['"]$/g, '');
        const rawValue = content.slice(separator + 1).trim();
        while (stack.length > 1 && indent <= stack[stack.length - 1].indent) stack.pop();
        const parent = stack[stack.length - 1].value;
        if (!parent || typeof parent !== 'object' || Array.isArray(parent)) throw new Error(`Invalid YAML nesting on line ${lineIndex + 1}.`);
        if (Object.prototype.hasOwnProperty.call(parent, key)) throw new Error(`Duplicate YAML key: ${key}`);
        if (rawValue === '') {
            parent[key] = {};
            stack.push({ indent, value: parent[key] });
        } else {
            parent[key] = parseYamlScalar(rawValue);
        }
    });
    return root;
}

function parseSopDocument(text, filename = '') {
    const trimmed = text.trim();
    if (!trimmed) throw new Error('SOP file is empty.');
    if (/\.json$/i.test(filename) || trimmed.startsWith('{')) return JSON.parse(trimmed);
    return parseSopYaml(trimmed);
}

function yamlScalar(value) {
    if (value === null || value === undefined) return 'null';
    if (typeof value === 'boolean' || typeof value === 'number') return String(value);
    if (typeof value === 'string' && /^[A-Za-z0-9_$.-]+$/.test(value)) return value;
    return JSON.stringify(String(value));
}

function serializeSopYaml(value, indent = 0) {
    return Object.entries(value).map(([key, item]) => {
        const prefix = `${' '.repeat(indent)}${key}:`;
        if (Array.isArray(item)) return `${prefix} [${item.map(yamlScalar).join(', ')}]`;
        if (item && typeof item === 'object') return `${prefix}\n${serializeSopYaml(item, indent + 2)}`;
        return `${prefix} ${yamlScalar(item)}`;
    }).join('\n');
}

function validateSopProfile(profile) {
    if (!profile || typeof profile !== 'object' || Array.isArray(profile)) throw new Error('SOP must be a JSON object.');
    if (profile.$schema !== SOP_SCHEMA) throw new Error(`SOP schema must be ${SOP_SCHEMA}.`);
    ['profile_id', 'sop_name', 'version', 'status'].forEach(key => {
        if (typeof profile[key] !== 'string' || !profile[key].trim()) throw new Error(`SOP field ${key} is required.`);
    });
    if (String(profile.unknown_rule_policy).toUpperCase() !== 'ERROR') throw new Error('unknown_rule_policy must be ERROR.');
    if (!profile.rules || typeof profile.rules !== 'object' || Array.isArray(profile.rules)) throw new Error('SOP rules must be an object.');
    const allowedEnforcement = new Set(['HARD_FAIL', 'WARNING', 'IGNORE', 'INFORMATIONAL']);
    Object.entries(profile.rules).forEach(([ruleId, level]) => {
        if (!ruleId || !allowedEnforcement.has(String(level).toUpperCase())) throw new Error(`Invalid SOP rule policy: ${ruleId}`);
    });
    const workflow = profile.workflow;
    if (!workflow || typeof workflow !== 'object') throw new Error('SOP workflow is required.');
    const validMethods = new Set(['feasibility_best', 'dp_v2_1_1', 'high_cai', 'gc_target', 'assembly_friendly']);
    if (!validMethods.has(workflow.design_method)) throw new Error('Unsupported SOP design_method.');
    if (workflow.design_mode !== undefined && !['baseline_1dp', 'paired_3dp'].includes(workflow.design_mode)) throw new Error('Unsupported SOP design_mode.');
    if (workflow.structure_settings !== undefined) {
        const settings = workflow.structure_settings;
        if (!settings || typeof settings !== 'object' || Array.isArray(settings)) throw new Error('Invalid structure_settings.');
        if (typeof settings.chain !== 'string' || !settings.chain.trim() || settings.chain.length > 8 || !Number.isInteger(settings.model) || settings.model < 1) throw new Error('Invalid structure chain/model settings.');
        if (!Array.isArray(settings.positions) || settings.positions.some(p => !Number.isSafeInteger(p) || p < 1) || new Set(settings.positions).size !== settings.positions.length) throw new Error('Invalid structure positions.');
        if (!Number.isFinite(settings.target_adaptiveness) || settings.target_adaptiveness <= 0 || settings.target_adaptiveness > 1 || !Number.isFinite(settings.penalty_weight) || settings.penalty_weight < 0 || settings.penalty_weight > 100) throw new Error('Invalid structure target/penalty settings.');
    }
    const requirements = workflow.sequence_requirements;
    if (!requirements || typeof requirements !== 'object') throw new Error('SOP sequence_requirements are required.');
    if (!Array.isArray(requirements.type_iis_enzymes) || requirements.type_iis_enzymes.some(name => !TYPE_IIS_PRESETS[name])) throw new Error('SOP contains an unsupported Type IIS enzyme.');
    if (!Array.isArray(requirements.custom_restriction_sites)) throw new Error('SOP custom_restriction_sites must be an array.');
    const policy = workflow.review_policy;
    if (!policy || typeof policy !== 'object') throw new Error('SOP review_policy is required.');
    ['cai', 'overall_gc', 'local_gc', 'type_iis', 'repeats', 'homopolymers', 'forbidden_motifs'].forEach(key => {
        if (!SOP_MODES.includes(String(policy[key]).toLowerCase())) throw new Error(`Invalid review policy mode: ${key}`);
    });
    return cloneJson(profile);
}

function loadSopProfile() {
    if (!ENABLE_CUSTOM_SOP) return cloneJson(DEFAULT_SOP_PROFILE);
    try {
        localStorage.removeItem('factorforge_active_sop_v1'); // Purge legacy cache with obsolete naming
        const saved = localStorage.getItem(SOP_STORAGE_KEY);
        if (!saved) return cloneJson(DEFAULT_SOP_PROFILE);
        const parsed = JSON.parse(saved);
        if (parsed && typeof parsed === 'object') {
            if (typeof parsed.sop_name === 'string') {
                parsed.sop_name = parsed.sop_name
                    .replace(/PlantForm\s+Doug\s+Biologics/gi, 'Therapeutic Biologics')
                    .replace(/PlantForm/gi, 'Therapeutic')
                    .replace(/Doug/gi, 'Balanced');
            }
            if (parsed.profile_id === 'sop-preset-b-plantform-biologics') {
                return cloneJson(SOP_PRESETS.preset_b.profile);
            }
        }
        return validateSopProfile(parsed);
    } catch (_) {
        localStorage.removeItem(SOP_STORAGE_KEY);
        return cloneJson(DEFAULT_SOP_PROFILE);
    }
}

function loadVersionedHistory() {
    try {
        const raw = JSON.parse(localStorage.getItem('factorforge_history') || '[]');
        if (Array.isArray(raw)) return raw.map(item => ({ ...item, schemaVersion: item.schemaVersion || 1 }));
        if (raw && Array.isArray(raw.items)) {
            return raw.items.map(item => ({
                ...item,
                schemaVersion: item.schemaVersion || raw.schemaVersion || 1,
            }));
        }
    } catch (_) {
        localStorage.removeItem('factorforge_history');
    }
    return [];
}

// State Management
const state = {
    sequence: '',
    engineMode: 'profile',
    objective: 'dp_v2_1_1',
    host: 'nbenthamiana',
    watermark: false,
    codonPolicy: 'plantform_balanced_codon_v1',
    saveDb: false,
    alignmentPage: 0,
    useTemplate: false,
    kozak: false,
    dinuc: false,
    customRestrictionSites: [],
    selectedTypeIisEnzymes: [],
    reviewerDisposition: null,
    results: null,
    reportData: null,
    isOptimizing: false,
    history: loadVersionedHistory()
};
state.activeSop = loadSopProfile();

// DOM Elements
const elements = {
    fileUpload: document.getElementById('fileUpload'),
    sequenceInput: document.getElementById('sequenceInput'),
    sequencePreview: document.getElementById('sequencePreview'),
    previewContainer: document.getElementById('previewContainer'),
    validationWarning: document.getElementById('validationWarning'),
    clearBtn: document.getElementById('clearBtn'),
    optimizeBtn: document.getElementById('optimizeBtn'),
    btnText: document.getElementById('btnText'),
    loadingIndicator: document.getElementById('loadingIndicator'),
    validationStatus: document.getElementById('validationStatus'),
    engineVersionSelect: document.getElementById('engineVersionSelect'),
    engineVersionDescription: document.getElementById('engineVersionDescription'),
    engineBadge: document.getElementById('engineBadge'),
    toggleWatermark: document.getElementById('toggleWatermark'),
    codonPolicyRadios: document.getElementsByName('codonPolicy'),
    activePolicyBadge: document.getElementById('activePolicyBadge'),
    designWorkspace: document.getElementById('designWorkspace'),
    designBriefPanel: document.getElementById('designBriefPanel'),
    resultsPanel: document.getElementById('resultsPanel'),
    emptyState: document.getElementById('emptyState'),
    resultsContainer: document.getElementById('resultsContainer'),
    caiValue: document.getElementById('caiValue'),
    gcValue: document.getElementById('gcValue'),
    polyaValue: document.getElementById('polyaValue'),
    hostProfileValue: document.getElementById('hostProfileValue'),
    optimizedSequence: document.getElementById('optimizedSequence'),
    jsonDetails: document.getElementById('jsonDetails'),
    downloadEvidencePackageZip: document.getElementById('downloadEvidencePackageZip'),
    downloadFasta: document.getElementById('downloadFasta'),
    exportInteractiveReport: document.getElementById('exportInteractiveReport'),
    downloadGenbank: document.getElementById('downloadGenbank'),
    copyBtn: document.getElementById('copyBtn'),
    constructIdRow: document.getElementById('constructIdRow'),
    constructIdDisplay: document.getElementById('constructIdDisplay'),
    copyConstructId: document.getElementById('copyConstructId'),
    submitValidationBtn: document.getElementById('submitValidationBtn'),
    alphafoldLink: document.getElementById('alphafoldLink'),
    esmatlasFoldLink: document.getElementById('esmatlasFoldLink'),
    copyJsonBtn: document.getElementById('copyJsonBtn'),
    statTotalCodons: document.getElementById('statTotalCodons'),
    statShifts: document.getElementById('statShifts'),
    statUnchanged: document.getElementById('statUnchanged'),
    statSubstitutions: document.getElementById('statSubstitutions'),
    statShiftBreakdown: document.getElementById('statShiftBreakdown'),
    statAaIdentity: document.getElementById('statAaIdentity'),
    interactiveTrackContainer: document.getElementById('interactiveTrackContainer'),
    canvasOriginalTrack: document.getElementById('canvasOriginalTrack'),
    canvasOptimizedTrack: document.getElementById('canvasOptimizedTrack'),
    trackAvailability: document.getElementById('trackAvailability'),
    inspectorPosition: document.getElementById('inspectorPosition'),
    inspectorCodons: document.getElementById('inspectorCodons'),
    inspectorRationale: document.getElementById('inspectorRationale'),
    inspectorConstraints: document.getElementById('inspectorConstraints'),
    toggleDetails: document.getElementById('toggleDetails'),
    detailsContent: document.getElementById('detailsContent'),
    toggleArrow: document.getElementById('toggleArrow'),
    themeToggle: document.getElementById('themeToggle'),
    themeIcon: document.getElementById('themeIcon'),
    objectiveRadios: document.getElementsByName('objective'),
    engineModeRadios: document.getElementsByName('engineMode'),
    engineSelector: document.getElementById('engineSelector'),
    hostSelect: document.getElementById('hostSelect'),
    saveDbToggle: document.getElementById('saveDbToggle'),
    saveDbStatus: document.getElementById('saveDbStatus'),
    databaseCapability: document.getElementById('databaseCapability'),
    appliedPolicySummary: document.getElementById('appliedPolicySummary'),
    resultContextSummary: document.getElementById('resultContextSummary'),
    designContractSummary: document.getElementById('designContractSummary'),
    comparisonDashboard: document.getElementById('comparisonDashboard'),
    comparisonNotice: document.getElementById('comparisonNotice'),
    comparisonMatrixBody: document.getElementById('comparisonMatrixBody'),
    provenanceBadges: document.getElementById('provenanceBadges'),
    codonAlignmentViewer: document.getElementById('codonAlignmentViewer'),
    alignmentRange: document.getElementById('alignmentRange'),
    alignmentPrev: document.getElementById('alignmentPrev'),
    alignmentNext: document.getElementById('alignmentNext'),
    implementedObjectives: document.getElementById('implementedObjectives'),
    experimentalObjectives: document.getElementById('experimentalObjectives'),
    packagedReferenceAssets: document.getElementById('packagedReferenceAssets'),
    useTemplateCheck: document.getElementById('useTemplate'),
    kozakToggle: document.getElementById('toggleKozak'),
    dinucToggle: document.getElementById('toggleDinuc'),
    customRestrictionSites: document.getElementById('customRestrictionSites'),
    optimizationSeed: document.getElementById('optimizationSeed'),
    typeIisEnzymes: document.getElementsByName('typeIisEnzyme'),
    inputLenBadge: document.getElementById('inputLenBadge'),
    inputGCBadge: document.getElementById('inputGCBadge'),
    origLen: document.getElementById('origLen'),
    optLen: document.getElementById('optLen'),
    origGC: document.getElementById('origGC'),
    optGCComp: document.getElementById('optGCComp'),
    origCAI: document.getElementById('origCAI'),
    optCAIComp: document.getElementById('optCAIComp'),
    mutationRate: document.getElementById('mutationRate'),
    aaIdentity: document.getElementById('aaIdentity'),
    candidateComparisonContainer: document.getElementById('candidateComparisonContainer'),
    candidateComparisonBody: document.getElementById('candidateComparisonBody'),
    customRestrictionResults: document.getElementById('customRestrictionResults'),
    customRestrictionResultsBody: document.getElementById('customRestrictionResultsBody'),
    mfeWarningBanner: document.getElementById('mfeWarningBanner'),
    gcTargetRange: document.getElementById('gcTargetRange'),
    resultsReport: document.getElementById('resultsReport'),
    resultsReportBody: document.getElementById('resultsReportBody'),
    gcChart: document.getElementById('gcChart'),
    gcZoneLabel: document.getElementById('gcZoneLabel'),
    historyList: document.getElementById('historyList'),
    clearHistory: document.getElementById('clearHistory'),
    // Guided workspace elements
    designWorkspace: document.getElementById('designWorkspace'),
    btnEditSettings: document.getElementById('btnEditSettings'),
    reviewConfigSummary: document.getElementById('reviewConfigSummary'),
    tabBtnSummary: document.getElementById('tabBtnSummary'),
    tabBtnComparison: document.getElementById('tabBtnComparison'),
    tabBtnDifferences: document.getElementById('tabBtnDifferences'),
    tabBtnTrace: document.getElementById('tabBtnTrace'),
    tabPaneComparison: document.getElementById('tabPaneComparison'),
    tabPaneTrace: document.getElementById('tabPaneTrace'),
    stepInputBadge: document.getElementById('stepInputBadge'),
    stepPolicyBadge: document.getElementById('stepPolicyBadge'),
    stepResultsBadge: document.getElementById('stepResultsBadge'),
    toggleArgSerDistribution: document.querySelector('input[name="toggleArgSerDistribution"]'),
    toggleSpliceScan: document.querySelector('input[name="toggleSpliceScan"]'),
    toggleFivePrimeBoundary: document.querySelector('input[name="toggleFivePrimeBoundary"]'),

    changelogBtn: document.getElementById('changelogBtn'),
    changelogModal: document.getElementById('changelogModal'),
    closeModal: document.getElementById('closeModal'),
    modalOverlay: document.getElementById('modalOverlay'),
    linkoutConsentModal: document.getElementById('linkoutConsentModal'),
    linkoutConsentOverlay: document.getElementById('linkoutConsentOverlay'),
    linkoutConsentTitle: document.getElementById('linkoutConsentTitle'),
    linkoutConsentBody: document.getElementById('linkoutConsentBody'),
    linkoutConsentClose: document.getElementById('linkoutConsentClose'),
    linkoutConsentCancel: document.getElementById('linkoutConsentCancel'),
    linkoutConsentContinue: document.getElementById('linkoutConsentContinue'),
    inputTypeBadge: document.getElementById('inputTypeBadge'),
    toastContainer: document.getElementById('toastContainer'),
    logoIcon: document.getElementById('logoIcon'),
    logoTitle: document.getElementById('logoTitle'),
    criterionCaiMode: document.getElementById('criterionCaiMode'),
    criterionGcMode: document.getElementById('criterionGcMode'),
    criterionLocalGcMode: document.getElementById('criterionLocalGcMode'),
    criterionTypeIisMode: document.getElementById('criterionTypeIisMode'),
    criterionRepeatsMode: document.getElementById('criterionRepeatsMode'),
    criterionHomopolymerMode: document.getElementById('criterionHomopolymerMode'),
    criterionMotifsMode: document.getElementById('criterionMotifsMode'),
    sopPresetSelect: document.getElementById('sopPresetSelect'),
    sopPresetTag: document.getElementById('sopPresetTag'),
    sopPresetDescription: document.getElementById('sopPresetDescription'),
    downloadSopExample: document.getElementById('downloadSopExample'),
    sopProfileName: document.getElementById('sopProfileName'),
    sopProfileMeta: document.getElementById('sopProfileMeta'),
    sopProfileBadge: document.getElementById('sopProfileBadge'),
    sopUploadZone: document.getElementById('sopUploadZone'),
    downloadSopTemplate: document.getElementById('downloadSopTemplate'),
    uploadSopButton: document.getElementById('uploadSopButton'),
    resetSopProfile: document.getElementById('resetSopProfile'),
    sopFileUpload: document.getElementById('sopFileUpload'),
    automatedDecisionValue: document.getElementById('automatedDecisionValue'),
    automatedDecisionSummary: document.getElementById('automatedDecisionSummary'),
    qcDecisionMatrix: document.getElementById('qcDecisionMatrix'),
    qcDecisionMatrixBody: document.getElementById('qcDecisionMatrixBody'),
    reviewerDisposition: document.getElementById('reviewerDisposition'),
    reviewerReason: document.getElementById('reviewerReason'),
    saveReviewerDisposition: document.getElementById('saveReviewerDisposition'),
    reviewerDispositionStatus: document.getElementById('reviewerDispositionStatus')
};

let chartInstance = null;

// Analytics helpers
function seqLenBucket(len) {
    if (len < 100) return '<100';
    if (len < 300) return '100-300';
    if (len < 1000) return '300-1000';
    return '>1000';
}
function caiBucket(cai) {
    if (cai < 0.7) return '<0.7';
    if (cai < 0.8) return '0.7-0.8';
    if (cai < 0.9) return '0.8-0.9';
    return '>0.9';
}
function gcBucket(gc, hostId = state.host) {
    // Host-aware telemetry bucketing (v3.3.0) — boundaries follow
    // the resolved reference band for the host instead of a fixed 55-65
    // assumption, so BY-2 and N. benthamiana don't get mislabeled buckets.
    const { gc_min, gc_max } = getGcRange(hostId);
    if (gc < gc_min) return `<${gc_min}`;
    if (gc <= gc_max) return `${gc_min}-${gc_max}`;
    return `>${gc_max}`;
}
function trackEvent(name, data) {
    try { window.va?.('event', { name, data }); } catch (_) {}
}

// Initialization
document.addEventListener('DOMContentLoaded', async () => {
    initTheme();
    applyStaticLabelPatches();
    applySopToUi(state.activeSop, { persist: false });
    initEventListeners();
    setupDesignPanel();
    handleSequenceChange({ target: elements.sequenceInput });
    await loadApiMetadata();
    applySopToUi(state.activeSop, { persist: false });
    updateDesignBriefSummary();


    renderHistory();
    console.log('FactorForge ready');
});

// Move existing controls, rather than recreating them: values and event handlers
// stay intact, including custom profiles retained in this browser.
function setupDesignPanel() {
    // Preserve the legacy DOM/hooks; hide custom-policy management after setup.

    const advanced = document.querySelector('#manualSopOverrides > div');
    advanced.prepend(document.getElementById('designObjectivePolicy'));
    document.querySelector('#manualSopOverrides > summary').textContent = 'Advanced settings';
    const experimental = document.createElement('details');
    experimental.id = 'experimentalSettings';
    experimental.className = 'disclosure-card mb-6 rounded-xl border border-slate-200 dark:border-slate-800 p-4';
    const summary = document.createElement('summary');
    summary.textContent = 'Experimental features';
    summary.className = 'cursor-pointer text-xs font-bold';
    experimental.append(summary);
    experimental.append(document.getElementById('engineSelector'));
    experimental.append(document.getElementById('experimentalObjectives'));
    experimental.append(document.getElementById('toggleWatermark').closest('label'));
    advanced.append(document.getElementById('useTemplate').closest('label'));
    document.getElementById('manualSopOverrides').after(experimental);
    const exampleLabels = ['FactorForge Default (example)', 'B · Protein design example',
        'C · Modular cloning example', 'D · Expression workflow example', 'E · Synthesis review example'];
    Array.from(elements.sopPresetSelect.options).forEach((option, index) => {
        if (exampleLabels[index]) option.textContent = exampleLabels[index];
    });
    const customOption = document.createElement('option');
    customOption.value = 'custom';
    customOption.textContent = 'Uploaded custom profile';
    customOption.disabled = true;
    elements.sopPresetSelect.append(customOption);
    updatePresetUiFromProfile(state.activeSop);
    const card = document.getElementById('sopProfileCard');
    const management = document.createElement('details');
    management.id = 'profileManagement';
    management.className = 'disclosure-card text-xs';
    const managementSummary = document.createElement('summary');
    managementSummary.className = 'cursor-pointer font-bold py-2';
    managementSummary.textContent = 'Upload, download or edit profile';
    management.append(managementSummary);
    management.append(elements.uploadSopButton.parentElement);
    management.append(document.getElementById('defaultSopSummary'));
    management.append(card.querySelector('a[target="_blank"]').parentElement);
    card.append(management);
    document.querySelectorAll('#designBriefPanel .info-icon').forEach(icon => {
        icon.tabIndex = 0;
        icon.setAttribute('aria-label', icon.nextElementSibling?.textContent.trim() || 'Setting details');
    });
    elements.uploadSopButton.parentElement.classList.add('flex-wrap');
    elements.downloadSopExample.title = 'Download an editable example configuration, not a validated laboratory SOP.';
    document.getElementById('sopProfileHeading').textContent = 'Sequence Design Profile';
    const note = document.createElement('p');
    note.id = 'profileExampleNote';
    note.className = 'text-xs text-slate-500 dark:text-slate-400';
    note.textContent = 'Example configuration · Adapt to your workflow';
    const boundary = document.createElement('p');
    boundary.className = 'mt-2 text-xs text-slate-500 dark:text-slate-400';
    boundary.textContent = 'These are example settings, not an approved laboratory SOP. Computational checks do not establish biological performance.';
    document.getElementById('defaultSopSummary').append(boundary);
    card.insertBefore(note, card.children[1]);
    elements.sopPresetSelect.setAttribute('aria-label', 'Design policy example');
    const resultNotice = document.createElement('p');
    resultNotice.id = 'settingsChangeNotice';
    resultNotice.className = 'hidden text-xs text-amber-700 dark:text-amber-300 mb-4';
    resultNotice.textContent = 'Settings changed. Existing results retain their original settings; generate again to apply the current profile.';
    document.getElementById('optimizeBtn').before(resultNotice);
    if (!ENABLE_CUSTOM_SOP) {
        management.hidden = true;
        elements.sopPresetSelect.parentElement.hidden = true;
        document.getElementById('manualSopOverrides').hidden = true;
        elements.sopFileUpload.disabled = true;
        elements.uploadSopButton.disabled = true;
        elements.downloadSopTemplate.disabled = true;
        document.getElementById('buildSopButton').disabled = true;
        elements.sopProfileName.textContent = 'FactorForge default policy';
        elements.sopProfileBadge.textContent = 'Default';
        note.textContent = 'Built-in FactorForge policy is applied automatically. Custom SOP upload, export and editing are paused.';
        document.getElementById('sopProfileHeading').textContent = 'Built-in sequence policy';
    }

}

// Loads server-owned GC ranges and validation labels. Supported hosts remain
// available to API/CLI clients, while the public web workflow stays focused on
// the production N. benthamiana path.
async function loadApiMetadata() {
    try {
        const response = await fetch(API_ENDPOINT, { method: 'GET', signal: AbortSignal.timeout(8000) });
        if (response.ok) {
            const data = await response.json();
            if (data.host_metadata && typeof data.host_metadata === 'object') {
                hostMetadata = data.host_metadata;
                Object.entries(data.host_metadata).forEach(([id, meta]) => {
                    if (meta && meta.gc_range) {
                        hostGcRanges[id] = meta.gc_range;
                    }
                });
            }
            if (Array.isArray(data.validation_checks)) {
                validationRegistry = data.validation_checks;
            }
            apiCapabilities = data.capabilities || {};
            const dpV21Available = Array.isArray(data.supported_objectives)
                && data.supported_objectives.includes('dp_v2_1_1');
            const dpV21Radio = document.getElementById('dpV21Radio');
            const dpV21Capability = document.getElementById('dpV21Capability');
            if (dpV21Radio) dpV21Radio.disabled = !dpV21Available;
            if (dpV21Capability) {
                dpV21Capability.textContent = dpV21Available
                    ? 'Available · engine 2.1.1-dev · calibration complete / holdout pending'
                    : 'Unavailable on this deployment';
            }
            const mlAvailable = Boolean(apiCapabilities.ml_preview?.available);
            elements.engineModeRadios.forEach(radio => {
                if (radio.value !== 'profile') radio.disabled = !mlAvailable;
            });
            elements.engineSelector?.classList.toggle('hidden', !mlAvailable);
            const dbAvailable = Boolean(apiCapabilities.db_save?.available);
            elements.saveDbToggle.disabled = !dbAvailable;
            elements.databaseCapability?.classList.toggle('hidden', !dbAvailable);
            elements.databaseCapability?.classList.toggle('flex', dbAvailable);
            elements.saveDbStatus.textContent = dbAvailable
                ? 'Enabled for this deployment'
                : 'Unavailable on this deployment';
        }
    } catch (_) {
        // Offline/dev fallback — getGcRange() uses OFFLINE_GC_RANGES;
        // validationRegistry stays [] instead of guessing server labels.
    }
}

const ENGINE_DESCRIPTIONS = {
    dp_v2_1_1: {
        badge: 'v2.1 recommended',
        badgeClass: 'bg-teal-50 text-teal-700 border-teal-100 dark:bg-teal-900/20 dark:text-teal-300 dark:border-teal-800/50',
        cardClass: 'bg-teal-50/80 dark:bg-teal-900/20 border-teal-100 dark:border-teal-800/40 text-teal-900 dark:text-teal-200',
        text: '<span class="font-bold text-teal-800 dark:text-teal-300">✨ Recommended:</span> 3-axis DP with position-dependent 5′ initiation ramp awareness (nt 1–45 clamped to 20–29% GC), homopolymer veto (&le;5 nt), global GC window, and host CAI maximization.'
    },
    feasibility_best: {
        badge: 'v2.0 standard',
        badgeClass: 'bg-emerald-50 text-emerald-700 border-emerald-100 dark:bg-emerald-900/20 dark:text-emerald-300 dark:border-emerald-800/50',
        cardClass: 'bg-emerald-50/80 dark:bg-emerald-900/20 border-emerald-100 dark:border-emerald-800/40 text-emerald-900 dark:text-emerald-200',
        text: '<span class="font-bold text-emerald-800 dark:text-emerald-300">Global Feasibility DP:</span> Classical exact dynamic programming enforcing host GC equilibrium and CAI maximization across the complete CDS.'
    },
    high_cai: {
        badge: 'rule v1.0',
        badgeClass: 'bg-blue-50 text-blue-700 border-blue-100 dark:bg-blue-900/20 dark:text-blue-300 dark:border-blue-800/50',
        cardClass: 'bg-blue-50/80 dark:bg-blue-900/20 border-blue-100 dark:border-blue-800/40 text-blue-900 dark:text-blue-200',
        text: '<span class="font-bold text-blue-800 dark:text-blue-300">Rule Profile:</span> Greedy selection favoring the highest codon adaptation index (CAI) for each residue.'
    },
    gc_target: {
        badge: 'rule v1.0',
        badgeClass: 'bg-slate-100 text-slate-700 border-slate-200 dark:bg-slate-800 dark:text-slate-300 dark:border-slate-700',
        cardClass: 'bg-slate-50 dark:bg-slate-900 border-slate-200 dark:border-slate-800 text-slate-800 dark:text-slate-200',
        text: '<span class="font-bold text-slate-800 dark:text-slate-300">Rule Profile:</span> Prioritizes overall CDS proximity to the host natural GC band.'
    },
    assembly_friendly: {
        badge: 'rule v1.0',
        badgeClass: 'bg-amber-50 text-amber-700 border-amber-100 dark:bg-amber-900/20 dark:text-amber-300 dark:border-amber-800/50',
        cardClass: 'bg-amber-50/80 dark:bg-amber-900/20 border-amber-100 dark:border-amber-800/40 text-amber-900 dark:text-amber-200',
        text: '<span class="font-bold text-amber-800 dark:text-amber-300">Rule Domesticator:</span> Automatically removes BsaI, BpiI, and BsmBI sites via synonymous substitution.'
    }
};

function updateEngineDescription(objective) {
    const meta = ENGINE_DESCRIPTIONS[objective] || ENGINE_DESCRIPTIONS.dp_v2_1_1;
    if (elements.engineBadge) {
        elements.engineBadge.textContent = meta.badge;
        elements.engineBadge.className = `shrink-0 px-2 py-1 rounded-lg border text-[9px] font-extrabold uppercase tracking-wider ${meta.badgeClass}`;
    }
    if (elements.engineVersionDescription) {
        elements.engineVersionDescription.innerHTML = meta.text;
        elements.engineVersionDescription.className = `p-3 border rounded-xl text-xs leading-relaxed ${meta.cardClass}`;
    }
}

function formatCustomSitesForEditor(sites) {
    return (sites || []).map(site => `${site.name}:${site.sequence}`).join('\n');
}

function updatePresetUiFromProfile(profile) {
    if (!elements.sopPresetSelect) return;
    const matchedKey = Object.keys(SOP_PRESETS).find(k => SOP_PRESETS[k].profile.profile_id === profile?.profile_id);
    if (matchedKey) {
        elements.sopPresetSelect.value = matchedKey;
        const preset = SOP_PRESETS[matchedKey];
        if (elements.sopPresetTag) {
            elements.sopPresetTag.textContent = matchedKey === 'preset_a' ? 'Default example active' : preset.tag;
            elements.sopPresetTag.className = 'px-2 py-0.5 rounded text-[9px] font-extrabold bg-emerald-100 text-emerald-800 border border-emerald-200 dark:bg-emerald-900/40 dark:text-emerald-300 dark:border-emerald-700';
        }
        if (elements.sopPresetDescription) {
            const descriptions = {
                preset_a: 'Default example: initiation GC constraints, homopolymer limits and configured Type IIS avoidance.',
                preset_b: 'Example with balanced codon distribution, Kozak settings and TpA suppression.',
                preset_c: 'Example with configured Type IIS constraints and MoClo assembly settings.',
                preset_d: 'Example with CAI-focused codon selection, Kozak settings and TpA suppression. Yield is not predicted.',
                preset_e: 'Example with codon-frequency distribution and synthesis-related checks. Vendor acceptance is not guaranteed.'
            };
            elements.sopPresetDescription.textContent = descriptions[matchedKey];
        }
        if (elements.downloadSopExample) {
            elements.downloadSopExample.href = preset.file;
            elements.downloadSopExample.download = '';
        }
    } else {
        elements.sopPresetSelect.value = 'custom';
        if (elements.sopPresetTag) {
            elements.sopPresetTag.textContent = 'Custom Policy';
            elements.sopPresetTag.className = 'px-2 py-0.5 rounded text-[9px] font-extrabold bg-amber-100 text-amber-800 border border-amber-200 dark:bg-amber-900/40 dark:text-amber-300 dark:border-amber-700';
        }
        if (elements.sopPresetDescription) {
            elements.sopPresetDescription.textContent = 'Custom sequence policy active. Click "Restore default" to reset to Preset A.';
        }
        if (elements.downloadSopExample) {
            elements.downloadSopExample.href = SOP_PRESETS.preset_a.file;
            elements.downloadSopExample.download = 'factorforge-sop-preset-a-conservative.yaml';
        }
    }
}

function applySopToUi(profile, { persist = true } = {}) {
    const validated = validateSopProfile(profile);
    validated.workflow.watermark = false;
    const workflow = validated.workflow;
    const requirements = workflow.sequence_requirements;
    const policy = workflow.review_policy;
    let objective = Array.from(elements.objectiveRadios).find(radio => radio.value === workflow.design_method && !radio.disabled);
    if (!objective) {
        objective = Array.from(elements.objectiveRadios).find(radio => !radio.disabled) || elements.objectiveRadios[0];
    }
    if (objective) {
        objective.checked = true;
        state.objective = objective.value;
        if (elements.engineVersionSelect) {
            elements.engineVersionSelect.value = objective.value;
            updateEngineDescription(objective.value);
        }
    }
    if (workflow.codon_policy) {
        state.codonPolicy = workflow.codon_policy;
        Array.from(elements.codonPolicyRadios).forEach(r => {
            r.checked = (r.value === workflow.codon_policy);
        });
        if (elements.activePolicyBadge) {
            elements.activePolicyBadge.textContent = workflow.codon_policy === 'plantform_balanced_codon_v1'
                ? 'Balanced Non-CpG'
                : (workflow.codon_policy === 'host_frequency' ? 'Host Frequency' : 'Max-CAI');
        }
    }
    if (workflow.watermark !== undefined && elements.toggleWatermark) {
        elements.toggleWatermark.checked = false;
        elements.toggleWatermark.disabled = true;
        state.watermark = false;
    }
    elements.optimizationSeed.value = requirements.reproducibility_seed ?? '';
    elements.kozakToggle.checked = Boolean(requirements.kozak_optimization);
    elements.dinucToggle.checked = Boolean(requirements.suppress_tpa);
    elements.useTemplateCheck.checked = Boolean(requirements.moclo_template);
    state.kozak = elements.kozakToggle.checked;
    state.dinuc = elements.dinucToggle.checked;
    state.useTemplate = elements.useTemplateCheck.checked;
    const selectedEnzymes = new Set(requirements.type_iis_enzymes);
    elements.typeIisEnzymes.forEach(input => { input.checked = selectedEnzymes.has(input.value); });
    elements.customRestrictionSites.value = formatCustomSitesForEditor(requirements.custom_restriction_sites);
    const modeMap = {
        cai: elements.criterionCaiMode, overall_gc: elements.criterionGcMode,
        local_gc: elements.criterionLocalGcMode, type_iis: elements.criterionTypeIisMode,
        repeats: elements.criterionRepeatsMode, homopolymers: elements.criterionHomopolymerMode,
        forbidden_motifs: elements.criterionMotifsMode
    };
    Object.entries(modeMap).forEach(([key, select]) => { select.value = String(policy[key]).toLowerCase(); });
    state.activeSop = validated;
    if (persist) localStorage.setItem(SOP_STORAGE_KEY, JSON.stringify(validated));
    const displayNames = {
        preset_a_conservative_baseline: 'FactorForge Default',
        preset_b_therapeutic_biologics: 'B · Protein design example',
        preset_c_moclo_goldenbraid: 'C · Modular cloning example',
        preset_d_high_yield_agro: 'D · Expression workflow example',
        preset_e_synthesis_feasibility: 'E · Synthesis review example'
    };
    elements.sopProfileName.textContent = displayNames[validated.profile_id] || validated.sop_name;
    elements.sopProfileMeta.textContent = `${validated.status === 'STANDARD_TEMPLATE' ? 'Default template' : (validated.status === 'PRESET_TEMPLATE' ? 'SOP Preset' : 'Custom SOP')} · v${validated.version} · stored only in this browser`;
    elements.sopProfileBadge.textContent = validated.status === 'STANDARD_TEMPLATE' ? 'Active' : (validated.status === 'PRESET_TEMPLATE' ? 'Preset' : 'Custom');
    if (!ENABLE_CUSTOM_SOP) {
        elements.sopProfileName.textContent = 'FactorForge default policy';
        elements.sopProfileBadge.textContent = 'Default';
    }
    updatePresetUiFromProfile(validated);
    if (!ENABLE_CUSTOM_SOP) {
        elements.sopProfileName.textContent = 'FactorForge default policy';
        elements.sopProfileBadge.textContent = 'Default';
    }
    const designMode = workflow.design_mode || 'baseline_1dp';
    state.researchDesignMode = designMode;
    document.querySelectorAll('input[name="researchDesignMode"]').forEach(input => { input.checked = input.value === designMode; });
    const structured = designMode === 'paired_3dp';
    document.getElementById('structureWebInputs').classList.toggle('hidden', !structured);
    document.getElementById('structureWebReviewed').checked = false;
    document.getElementById('structureWebFile').value = '';
    const settings = workflow.structure_settings || {chain:'A',model:1,positions:[],target_adaptiveness:0.4,penalty_weight:2.0};
    document.getElementById('structureWebChain').value = settings.chain;
    document.getElementById('structureWebModel').value = settings.model;
    document.getElementById('structureWebPositions').value = settings.positions.join(',');
    document.getElementById('structureWebTarget').value = settings.target_adaptiveness;
    document.getElementById('structureWebPenalty').value = settings.penalty_weight;
    elements.btnText.textContent = structured ? 'Generate 1DP / 3DP comparison' : 'Generate';
    document.getElementById('selectedWorkflowNote').textContent = structured ? '3DP settings loaded. Upload the matching structure and review positions before running.' : '1DP selected: sequence-only design in this form.';
    updateDesignBriefSummary();
    handleSequenceChange({target:elements.sequenceInput});
}

function captureSopFromUi() {
    if (!ENABLE_CUSTOM_SOP) {
        const profile = cloneJson(DEFAULT_SOP_PROFILE);
        profile.workflow.watermark = false;
        profile.workflow.design_method = state.researchDesignMode === 'paired_3dp' ? 'feasibility_best' : state.objective;
        if (state.researchDesignMode === 'paired_3dp') profile.workflow.design_mode = 'paired_3dp';
        state.activeSop = profile;
        return;
    }

    const profile = cloneJson(state.activeSop || DEFAULT_SOP_PROFILE);
    profile.status = profile.profile_id === DEFAULT_SOP_PROFILE.profile_id ? 'CUSTOMIZED_LOCAL' : profile.status;
    profile.workflow.design_method = state.objective;
    profile.workflow.codon_policy = state.codonPolicy;
    profile.workflow.watermark = false;
    if (state.researchDesignMode === 'paired_3dp') {
        profile.workflow.design_mode = 'paired_3dp';
        const raw = document.getElementById('structureWebPositions').value.trim();
        if (raw && raw.split(',').some(p => !/^\d+$/.test(p.trim()))) throw new Error('Structure positions must be comma-separated positive integers');
        profile.workflow.structure_settings = {chain:document.getElementById('structureWebChain').value.trim(),model:Number(document.getElementById('structureWebModel').value),positions:raw ? raw.split(',').map(Number) : [],target_adaptiveness:Number(document.getElementById('structureWebTarget').value),penalty_weight:Number(document.getElementById('structureWebPenalty').value)};
    } else {
        delete profile.workflow.design_mode;
        delete profile.workflow.structure_settings;
    }

    profile.workflow.sequence_requirements = {
        reproducibility_seed: elements.optimizationSeed.value.trim() === '' ? null : Number(elements.optimizationSeed.value),
        kozak_optimization: elements.kozakToggle.checked,
        suppress_tpa: elements.dinucToggle.checked,
        type_iis_enzymes: Array.from(elements.typeIisEnzymes).filter(input => input.checked).map(input => input.value),
        custom_restriction_sites: parseCustomRestrictionSites(elements.customRestrictionSites.value),
        moclo_template: elements.useTemplateCheck.checked
    };
    profile.workflow.review_policy = {
        cai: elements.criterionCaiMode.value.toUpperCase(), overall_gc: elements.criterionGcMode.value.toUpperCase(),
        local_gc: elements.criterionLocalGcMode.value.toUpperCase(), type_iis: elements.criterionTypeIisMode.value.toUpperCase(),
        repeats: elements.criterionRepeatsMode.value.toUpperCase(), homopolymers: elements.criterionHomopolymerMode.value.toUpperCase(),
        forbidden_motifs: elements.criterionMotifsMode.value.toUpperCase()
    };
    state.activeSop = validateSopProfile(profile);
    localStorage.setItem(SOP_STORAGE_KEY, JSON.stringify(state.activeSop));
    elements.sopProfileMeta.textContent = `Custom local settings · v${state.activeSop.version} · stored only in this browser`;
    elements.sopProfileBadge.textContent = 'Custom';
    updatePresetUiFromProfile(state.activeSop);
}

function downloadActiveSop() {
    if (!ENABLE_CUSTOM_SOP) return;
    captureSopFromUi();
    const safeId = state.activeSop.profile_id.replace(/[^a-z0-9_-]+/gi, '_');
    const exported = cloneJson(state.activeSop);
    exported.workflow.design_mode = state.researchDesignMode || 'baseline_1dp';
    const help = '# Selected settings snapshot. Edit values or change options in the UI, then upload again.\n# design_mode: baseline_1dp | paired_3dp\n# design_method (1DP): feasibility_best | dp_v2_1_1 | high_cai | gc_target | assembly_friendly\n# rules: HARD_FAIL (block) | WARNING (report) | INFORMATIONAL | IGNORE\n# review_policy: REQUIRED | PREFERRED | IGNORE\n# booleans: true | false; watermark is currently forced false.\n# 3DP structure file and review confirmation are never saved here; upload/review again.\n';
    downloadTextArtifact(`${help}${serializeSopYaml(exported)}\n`, 'application/yaml', `FactorForge_SOP_${safeId}_v${state.activeSop.version}.yaml`);
}

function applySopFile(file, input = null) {
    if (!ENABLE_CUSTOM_SOP) return;
    if (!file) return;
    if (file.size > 128 * 1024) {
        showToast('SOP file must be 128 KB or smaller.', 'error');
        if (input) input.value = '';
        return;
    }
    const reader = new FileReader();
    reader.onload = () => {
        try {
            applySopToUi(parseSopDocument(String(reader.result), file.name));
            showToast('SOP validated, applied, and saved in this browser.', 'success');
        } catch (error) {
            showToast(`SOP rejected: ${error.message}`, 'error');
        } finally {
            if (input) input.value = '';
        }
    };
    reader.onerror = () => showToast('Unable to read SOP file.', 'error');
    reader.readAsText(file);
}

function handleSopUpload(event) {
    applySopFile(event.target.files?.[0], event.target);
}

function restoreDefaultSop() {
    localStorage.removeItem(SOP_STORAGE_KEY);
    applySopToUi(cloneJson(DEFAULT_SOP_PROFILE), { persist: false });
    showToast('Default conservative SOP restored.', 'success');
}

function initEventListeners() {
    // Input Handling
    elements.fileUpload.addEventListener('change', handleFileUpload);
    elements.sequenceInput.addEventListener('input', debounce(handleSequenceChange, 300));
    elements.clearBtn.addEventListener('click', clearAll);
    elements.downloadSopTemplate.addEventListener('click', downloadActiveSop);
    elements.uploadSopButton.addEventListener('click', () => elements.sopFileUpload.click());
    elements.sopUploadZone.addEventListener('click', () => elements.sopFileUpload.click());
    ['dragenter', 'dragover'].forEach(type => elements.sopUploadZone.addEventListener(type, event => {
        event.preventDefault();
        elements.sopUploadZone.classList.add('ring-2', 'ring-emerald-500');
    }));
    ['dragleave', 'drop'].forEach(type => elements.sopUploadZone.addEventListener(type, event => {
        event.preventDefault();
        elements.sopUploadZone.classList.remove('ring-2', 'ring-emerald-500');
    }));
    elements.sopUploadZone.addEventListener('drop', event => applySopFile(event.dataTransfer?.files?.[0]));
    elements.sopFileUpload.addEventListener('change', handleSopUpload);
    elements.resetSopProfile.addEventListener('click', restoreDefaultSop);

    // Multi-SOP Presets Selector
    if (elements.sopPresetSelect) {
        elements.sopPresetSelect.addEventListener('change', (e) => {
            if (!ENABLE_CUSTOM_SOP) return;
            const presetKey = e.target.value;
            const targetPreset = SOP_PRESETS[presetKey];
            if (targetPreset) {
                applySopToUi(cloneJson(targetPreset.profile), { persist: true });
                showToast(`Applied ${e.target.selectedOptions[0].textContent}`, 'info');
            }
        });
    }

    // Engine / Version Selector
    if (elements.engineVersionSelect) {
        elements.engineVersionSelect.addEventListener('change', (e) => {
            state.objective = e.target.value;
            Array.from(elements.objectiveRadios).forEach(radio => {
                radio.checked = (radio.value === e.target.value);
            });
            updateEngineDescription(e.target.value);
            captureSopFromUi();
            updateDesignBriefSummary();
        });
    }

    elements.objectiveRadios.forEach(radio => {
        radio.addEventListener('change', (e) => {
            state.objective = e.target.value;
            if (elements.engineVersionSelect) {
                elements.engineVersionSelect.value = e.target.value;
                updateEngineDescription(e.target.value);
            }
            captureSopFromUi();
            updateDesignBriefSummary();
        });
    });

    if (elements.toggleWatermark) {
        elements.toggleWatermark.addEventListener('change', (e) => {
            e.target.checked = false;
            state.watermark = false;
            captureSopFromUi();
            updateDesignBriefSummary();
        });
    }

    if (elements.codonPolicyRadios) {
        elements.codonPolicyRadios.forEach(radio => {
            radio.addEventListener('change', (e) => {
                state.codonPolicy = e.target.value;
                if (elements.activePolicyBadge) {
                    elements.activePolicyBadge.textContent = e.target.value === 'plantform_balanced_codon_v1'
                        ? 'Balanced Non-CpG'
                        : (e.target.value === 'host_frequency' ? 'Host Frequency' : 'Max-CAI');
                }
                captureSopFromUi();
                updateDesignBriefSummary();
            });
        });
    }
    elements.engineModeRadios.forEach(radio => {
        radio.addEventListener('change', (e) => {
            state.engineMode = e.target.value;
            state.alignmentPage = 0;
        });
    });
    elements.hostSelect.addEventListener('change', (e) => {
        state.host = e.target.value;
        updateDesignBriefSummary();


    });
    elements.saveDbToggle.addEventListener('change', (e) => {
        state.saveDb = e.target.checked;
    });
    elements.alignmentPrev.addEventListener('click', () => changeAlignmentPage(-1));
    elements.alignmentNext.addEventListener('click', () => changeAlignmentPage(1));

    elements.useTemplateCheck.addEventListener('change', (e) => {
        state.useTemplate = e.target.checked;
        captureSopFromUi();
        updateDesignBriefSummary();


    });

    elements.kozakToggle.addEventListener('change', (e) => {
        state.kozak = e.target.checked;
        captureSopFromUi();
        updateDesignBriefSummary();


    });
    elements.dinucToggle.addEventListener('change', (e) => {
        state.dinuc = e.target.checked;
        captureSopFromUi();
        updateDesignBriefSummary();


    });
    elements.customRestrictionSites.addEventListener('input', () => {
        state.customRestrictionSites = [];
        updateDesignBriefSummary();


    });
    elements.customRestrictionSites.addEventListener('change', captureSopFromUi);
    elements.typeIisEnzymes.forEach(input => input.addEventListener('change', () => { captureSopFromUi(); updateDesignBriefSummary();

 }));
    document.querySelectorAll('.criterion-mode').forEach(select => select.addEventListener('change', captureSopFromUi));
    elements.optimizationSeed.addEventListener('change', captureSopFromUi);
    elements.saveReviewerDisposition.addEventListener('click', saveReviewerDisposition);
    elements.clearHistory.addEventListener('click', clearHistory);

    // Action
    elements.optimizeBtn.addEventListener('click', runOptimization);

    document.getElementById('structureWebFile')?.addEventListener('change', async (e) => {
        const file = e.target.files?.[0];
        if (!file) return;
        state.sampleStructureText = null;
        if (file.size > 250000) {
            document.getElementById('structureRecommendationStatus').textContent = 'File too large (> 250 KB). Please select a file under 250 KB.';
            return;
        }
        document.getElementById('structureRecommendationStatus').textContent = 'Analyzing structure coordinates and pause positions...';
        try {
            await recommendStructurePositions();
            if (structureSuggestion?.suggestion?.recommended_positions?.length) {
                document.getElementById('structureWebPositions').value = structureSuggestion.suggestion.recommended_positions.join(',');
                document.getElementById('structureWebReviewed').checked = true;
                document.getElementById('structureRecommendationStatus').textContent = `✅ Structure loaded. Recommended pause positions (${structureSuggestion.suggestion.recommended_positions.join(', ')}) auto-applied.`;
            }
        } catch (err) {
            console.warn('Auto recommend on upload error:', err);
        }
    });

    document.getElementById('recommendStructurePositions').addEventListener('click', recommendStructurePositions);
    document.getElementById('acceptStructurePositions').addEventListener('click',()=>{
        if (!structureSuggestion || structureSuggestion.revision!==structureRecommendationRevision) return;
        document.getElementById('structureWebPositions').value=structureSuggestion.suggestion.recommended_positions.join(',');
        document.getElementById('structureWebReviewed').checked=true;
        document.getElementById('structureRecommendationStatus').textContent='Suggested positions copied and ready for 3DP design.';
    });
    ['structureWebChain','structureWebModel','structureWebSource'].forEach(id=>document.getElementById(id).addEventListener('change',invalidateStructureSuggestion));
    elements.sequenceInput.addEventListener('input',invalidateStructureSuggestion);
    document.getElementById('structureWebPositions').addEventListener('input',()=>{document.getElementById('structureWebReviewed').checked=false;});

    const syncDesignModeUi = () => {
        const checkedRadio = document.querySelector('input[name="researchDesignMode"]:checked');
        if (!checkedRadio) return;
        const structured = checkedRadio.value === 'paired_3dp';
        state.researchDesignMode = checkedRadio.value;
        elements.btnText.textContent = structured ? 'Generate 3DP (Pilot)' : 'Generate 1DP';
        document.getElementById('structureWebInputs')?.classList.toggle('hidden', !structured);
        const noteEl = document.getElementById('selectedWorkflowNote');
        if (noteEl) {
            noteEl.textContent = structured
                ? '3DP (In Development): Experimental structure-conditioned design. Upload an exact-match PDB/mmCIF file.'
                : '1DP selected: standard production sequence design.';
        }
        updateDesignBriefSummary();
        handleSequenceChange({ target: elements.sequenceInput });
    };
    document.querySelectorAll('input[name="researchDesignMode"]').forEach(input => {
        input.addEventListener('change', syncDesignModeUi);
    });
    syncDesignModeUi();

    document.getElementById('openLocal3dViewerBtn')?.addEventListener('click', () => {
        const container = document.getElementById('structure3dViewerContainer');
        if (container) {
            container.classList.remove('hidden');
            container.scrollIntoView({ behavior: 'smooth', block: 'center' });
        }
    });

    // Results Actions
    elements.downloadEvidencePackageZip?.addEventListener('click', downloadEvidencePackageZip);
    elements.downloadFasta.addEventListener('click', () => downloadFile('fasta'));
    elements.exportInteractiveReport?.addEventListener('click', () => {
        if (!state.reportData) return showToast('Run an optimization before exporting a report.', 'info');
        downloadInteractiveReport(state.reportData);
    });
    elements.downloadGenbank.addEventListener('click', () => downloadFile('genbank'));
    elements.copyBtn.addEventListener('click', copyToClipboard);
    elements.copyConstructId.addEventListener('click', copyConstructId);
    elements.submitValidationBtn.addEventListener('click', submitValidation);
    elements.copyJsonBtn.addEventListener('click', copyJson);
    elements.toggleDetails.addEventListener('click', toggleDetailsPanel);
    elements.themeToggle.addEventListener('click', toggleTheme);
    elements.changelogBtn.addEventListener('click', toggleChangelog);
    initWorkspaceTabs();
    elements.closeModal.addEventListener('click', toggleChangelog);
    elements.modalOverlay.addEventListener('click', toggleChangelog);
    elements.logoIcon.addEventListener('click', reloadPage);
    elements.logoTitle.addEventListener('click', reloadPage);

    // Third-party structure linkout — gate navigation behind an explicit
    // consent step naming the actual external operator (no result yet → no-op).
    elements.alphafoldLink.addEventListener('click', (e) => {
        e.preventDefault();
        if (elements.alphafoldLink.getAttribute('href') === '#') return;
        openLinkoutConsent('alphafold', elements.alphafoldLink.href);
    });
    elements.esmatlasFoldLink.addEventListener('click', (e) => {
        e.preventDefault();
        if (elements.esmatlasFoldLink.getAttribute('href') === '#') return;
        openLinkoutConsent('esmatlas', elements.esmatlasFoldLink.href);
    });
    elements.linkoutConsentCancel.addEventListener('click', closeLinkoutConsent);
    elements.linkoutConsentClose.addEventListener('click', closeLinkoutConsent);
    elements.linkoutConsentOverlay.addEventListener('click', closeLinkoutConsent);
    elements.linkoutConsentContinue.addEventListener('click', closeLinkoutConsent);
}

function reloadPage() {
    window.location.reload();
}

function applyStaticLabelPatches() {
    if (!elements.submitValidationBtn) return;
    const label = elements.submitValidationBtn.children[1];
    if (label) label.textContent = 'Share Wet-lab Results (GitHub)';
}

function updateDesignBriefSummary() {
    if (!elements.appliedPolicySummary) return;
    const host = elements.hostSelect?.selectedOptions?.[0]?.textContent?.trim() || 'N. benthamiana';
    const selectedObjective = Array.from(elements.objectiveRadios).find(radio => radio.checked)?.value || state.objective;
    const methodLabels = {
        dp_v2_1_1: 'DP v2.1.1 (3-axis initiation)',
        feasibility_best: 'DP v2.0.1 (global feasibility)',
        high_cai: 'Rule · High CAI',
        gc_target: 'Rule · GC Target',
        assembly_friendly: 'Rule · Assembly Friendly'
    };
    const requirements = [];
    if (elements.useTemplateCheck?.checked) requirements.push('MoClo');
    if (state.watermark) requirements.push('watermark twin');
    if (state.codonPolicy === 'plantform_balanced_codon_v1') requirements.push('Balanced non-CpG');
    else if (state.codonPolicy === 'host_frequency') requirements.push('Host-frequency codons');
    else if (state.codonPolicy === 'max_cai') requirements.push('CAI-focused codons');
    if (elements.optimizationSeed?.value.trim() !== '') requirements.push(`seed ${elements.optimizationSeed.value.trim()}`);
    const enzymes = Array.from(elements.typeIisEnzymes || [])
        .filter(input => input.checked)
        .map(input => input.value);
    if (enzymes.length) requirements.push(`avoid ${enzymes.join('/')}`);
    if (elements.customRestrictionSites?.value.trim()) requirements.push('custom restriction sites');
    if (elements.kozakToggle?.checked) requirements.push('Kozak handling');
    if (elements.dinucToggle?.checked) requirements.push('TpA reduction');
    const requirementText = requirements.length ? ` · ${requirements.join(' · ')}` : '';
    const method = methodLabels[selectedObjective] || 'deterministic design';
    elements.appliedPolicySummary.textContent = state.researchDesignMode === 'paired_3dp'
        ? '3DP selected · direct experimental testing · uploaded structure and reviewed positions required · approval/evidence recording deferred'
        : `1DP selected · ${host} · ${method}${requirementText}`;
    const notice = document.getElementById('settingsChangeNotice');
    if (notice && state.results && !state.isOptimizing) notice.classList.remove('hidden');
}

function isProteinInputResult(res) {
    const explicitType = [
        res?.validation?.input_type,
        res?.input_type,
        res?.sequence_type
    ].find(Boolean);

    if (explicitType) {
        const normalized = String(explicitType).toLowerCase();
        return normalized === 'protein' || normalized === 'amino_acid' || normalized === 'amino-acid';
    }

    const seq = state.sequence || '';
    const isDNA = /^[ACGT]+$/i.test(seq);
    const isProtein = /^[ACDEFGHIKLMNPQRSTVWY*]+$/i.test(seq);
    return isProtein && !isDNA;
}

// Theme Handling
function initTheme() {
    const savedTheme = localStorage.getItem('theme') || 'light';
    if (savedTheme === 'dark') {
        document.documentElement.classList.add('dark');
        elements.themeIcon.textContent = '☀️';
    }
}

function toggleTheme() {
    const isDark = document.documentElement.classList.toggle('dark');
    localStorage.setItem('theme', isDark ? 'dark' : 'light');
    elements.themeIcon.textContent = isDark ? '☀️' : '🌙';
    showToast(`${isDark ? 'Dark' : 'Light'} mode enabled`, 'info');

    // Re-render chart if it exists to update colors
    if (state.results) renderGCGraph(state.results.optimized_sequence, getResultHostProfile(state.results));
}

// Handler Functions
function handleFileUpload(e) {
    const file = e.target.files[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (event) => {
        const content = event.target.result;
        elements.sequenceInput.value = content;
        handleSequenceChange({ target: { value: content } });
        showToast('File uploaded successfully', 'success');
    };
    reader.onerror = () => showToast('Error reading file', 'error');
    reader.readAsText(file);
}

function handleSequenceChange(e) {
    let rawValue = e.target.value;

    // Normalize exactly one FASTA record and preserve invalid characters so
    // they can be reported instead of silently deleted.
    const lines = rawValue.replace(/\r/g, '').split('\n').map(line => line.trim()).filter(Boolean);
    const headerIndexes = lines.map((line, index) => line.startsWith('>') ? index : -1).filter(index => index >= 0);
    let sequenceOnly = '';
    let parseError = '';
    if (headerIndexes.length > 1) {
        parseError = 'Multiple FASTA records detected. Upload one sequence at a time.';
    } else if (headerIndexes.length === 1 && headerIndexes[0] !== 0) {
        parseError = 'FASTA header must be the first non-empty line.';
    } else {
        sequenceOnly = (headerIndexes.length === 1 ? lines.slice(1) : lines).join('');
        sequenceOnly = sequenceOnly.replace(/\s/g, '').toUpperCase();
    }

    // Detection Regex
    const dnaRegex = /^[ACGTacgt]+$/;
    const proteinRegex = /^[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy*]+$/;

    const isDNA = dnaRegex.test(sequenceOnly);
    const isProtein = proteinRegex.test(sequenceOnly);

    // Update Badge & Warning
    const dnaInputWarning = document.getElementById('dnaInputWarning');
    if (sequenceOnly.length > 0 || parseError) {
        elements.inputTypeBadge.classList.remove('hidden');
        if (parseError) {
            elements.inputTypeBadge.textContent = 'Invalid FASTA';
            elements.inputTypeBadge.className = 'ml-2 px-1.5 py-0.5 rounded bg-rose-100 text-rose-700 font-bold text-[9px] uppercase';
            elements.validationWarning.classList.remove('hidden');
            elements.validationWarning.innerHTML = `<span class="mr-2">⚠️</span> ${parseError}`;
            if (dnaInputWarning) dnaInputWarning.classList.add('hidden');
        } else if (isDNA) {
            elements.inputTypeBadge.textContent = 'DNA';
            elements.inputTypeBadge.className = 'ml-2 px-1.5 py-0.5 rounded bg-blue-50 dark:bg-blue-900/20 text-blue-600 dark:text-blue-400 font-semibold text-[9px] uppercase';
            elements.validationWarning.classList.add('hidden');
            if (dnaInputWarning) dnaInputWarning.classList.remove('hidden');
        } else if (isProtein) {
            elements.inputTypeBadge.textContent = 'Protein';
            elements.inputTypeBadge.className = 'ml-2 px-1.5 py-0.5 rounded bg-emerald-50 dark:bg-emerald-900/20 text-emerald-600 dark:text-emerald-400 font-semibold text-[9px] uppercase';
            elements.validationWarning.classList.add('hidden');
            if (dnaInputWarning) dnaInputWarning.classList.add('hidden');
        } else {
            elements.inputTypeBadge.textContent = 'Mixed/Invalid';
            elements.inputTypeBadge.className = 'ml-2 px-1.5 py-0.5 rounded bg-rose-100 text-rose-700 font-bold text-[9px] uppercase';
            elements.validationWarning.classList.remove('hidden');
            const invalid = [...new Set(sequenceOnly.split('').filter(char => !/[ACGTNacgtn*ACDEFGHIKLMNPQRSTVWY]/.test(char)))].join(', ');
            elements.validationWarning.innerHTML = `<span class="mr-2">⚠️</span> Invalid characters: ${escapeHtml(invalid || 'input')}`;
            if (dnaInputWarning) dnaInputWarning.classList.add('hidden');
        }
    } else {
        elements.inputTypeBadge.classList.add('hidden');
        elements.validationWarning.classList.add('hidden');
        if (dnaInputWarning) dnaInputWarning.classList.add('hidden');
    }

    if (state.sequence !== (parseError ? '' : sequenceOnly)) invalidateStructureSuggestion();
    state.sequence = parseError ? '' : sequenceOnly;
    const validInput = !parseError && sequenceOnly.length >= 3 && (isDNA || (isProtein && !isDNA));
    elements.optimizeBtn.disabled = !validInput || state.isOptimizing;

    // Update Stats (only treat as protein if it's NOT also valid DNA)
    updateInputStats(state.sequence, isProtein && !isDNA);

    // Update Preview
    if (state.sequence.length > 0) {
        elements.previewContainer.classList.remove('hidden');
        elements.sequencePreview.textContent = state.sequence.substring(0, 150) + (state.sequence.length > 150 ? '...' : '');
    } else {
        elements.previewContainer.classList.add('hidden');
    }
}

function updateInputStats(seq, isProtein = false) {
    const len = seq.length;
    if (len > 0) {
        elements.inputLenBadge.textContent = `${len} ${isProtein ? 'aa' : 'bp'}`;
        elements.inputLenBadge.style.opacity = "1";

        if (isProtein) {
            elements.inputGCBadge.textContent = "GC: N/A";
            elements.inputGCBadge.style.opacity = "0.5";
            elements.inputGCBadge.title = "Not applicable for protein sequences";
        } else {
            const gc = calculateGC(seq);
            elements.inputGCBadge.textContent = `GC: ${gc}%`;
            elements.inputGCBadge.style.opacity = "1";
            elements.inputGCBadge.removeAttribute("title");
        }
    } else {
        elements.inputLenBadge.textContent = "0 bp";
        elements.inputLenBadge.style.opacity = "0.4";
        elements.inputGCBadge.textContent = "GC: 0%";
        elements.inputGCBadge.style.opacity = "0.4";
    }
}

let structureSuggestion = null;
let structureRecommendationRevision = 0;
function invalidateStructureSuggestion() {
    structureRecommendationRevision++;
    structureSuggestion = null;
    document.getElementById('acceptStructurePositions').disabled = true;
    document.getElementById('structureWebReviewed').checked = false;
    document.getElementById('structureRecommendationStatus').textContent = 'Input changed. Review new suggestions or enter positions manually.';
}
async function recommendStructurePositions() {
    let revision = structureRecommendationRevision;
    const button = document.getElementById('recommendStructurePositions');
    button.disabled = true;
    try {
        handleSequenceChange({target:elements.sequenceInput});
        revision = structureRecommendationRevision;
        const file = document.getElementById('structureWebFile')?.files?.[0] || null;
        const validFile = file && file.size <= 250000 && /\.(pdb|cif|mmcif)$/i.test(file.name);
        if (!validFile) throw new Error('Upload an exact-match PDB/mmCIF file before requesting positions.');
        if (!apiCapabilities.structure_recommendation?.available) throw new Error('Recommendations are unavailable on this backend');
        const source = document.getElementById('structureWebSource').value;
        const response = await fetch(API_ENDPOINT,{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({action:'recommend_structure',sequence:state.sequence,structure_context:{content:await file.text(),format:/\.pdb$/i.test(file.name)?'pdb':'mmcif',chain:document.getElementById('structureWebChain').value.trim(),model:Number(document.getElementById('structureWebModel').value),is_alphafold:source==='unknown'?null:source==='alphafold'}})});
        const data = await response.json();
        if (!response.ok || !data.recommendation) throw new Error(data.error || 'Recommendation failed');
        if (revision!==structureRecommendationRevision) return;
        const suggestion = data.recommendation;
        const ready = suggestion.status==='AVAILABLE' && suggestion.recommended_positions.length>0 && suggestion.reviewed===false;
        structureSuggestion = ready ? {revision,suggestion} : null;
        document.getElementById('acceptStructurePositions').disabled = !ready;
        document.getElementById('structureWebReviewed').checked = false;
        document.getElementById('structureRecommendationStatus').textContent = ready ? `Suggested protein positions: ${suggestion.recommended_positions.join(', ')}. Source: ${suggestion.provenance_basis}. Review before accepting; these are observation-based hypotheses, not confirmed pause/linker sites.` : `${suggestion.status}: ${suggestion.reason || 'No supported region'}. Manual reviewed positions remain available.`;
    } catch(error) {
        if (revision !== structureRecommendationRevision) return;
        structureSuggestion = null;
        document.getElementById('acceptStructurePositions').disabled = true;
        document.getElementById('structureWebReviewed').checked = false;
        document.getElementById('structureRecommendationStatus').textContent = error.message;
    } finally { button.disabled=false; }
}
function scientificReviewHtml(model) {
    const entries = model.candidate_scientific_reviews?.length ? model.candidate_scientific_reviews : [{label:'Selected candidate',scientific_review:model.scientific_review}];
    const rows = entries.filter(e=>e.scientific_review).map(entry=>{
        const observation=entry.scientific_review;
        const scan=observation.rna_risks || {};
        const mfe=observation.initiation_mfe || {};
        const value=mfe.status==='AVAILABLE' && Number.isFinite(mfe.mfe_kcal_mol)?`${mfe.mfe_kcal_mol} kcal/mol`:'Not measured';
        const count=value=>Number.isFinite(value)?String(value):'N/A';
        return `<div class="mt-2 rounded border border-slate-300 p-3 dark:border-slate-600"><strong>${escapeHtml(entry.label || entry.id)}</strong><p>RNA scan: ${escapeHtml(scan.summary_status || 'NOT_RECORDED')} · donor ${count(scan.donor_count)}, polyA ${count(scan.polya_count)}, ARE ${count(scan.are_count)}</p><p>Initiation MFE: ${escapeHtml(mfe.status || 'NOT_RECORDED')} · ${escapeHtml(value)} · ${escapeHtml(mfe.input_scope || 'unknown scope')} · ${escapeHtml(mfe.engine || 'unknown engine')}</p><p>Window [${escapeHtml(String(mfe.window_start ?? '?'))}, ${escapeHtml(String(mfe.window_end ?? '?'))}) · ${escapeHtml(String(mfe.temperature_c ?? '?'))} °C. ${escapeHtml(mfe.reason || mfe.advisory_flag || '')}</p></div>`;
    }).join('');
    return rows ? `<section class="mt-4 text-xs" aria-label="RNA risk and initiation observations"><h3 class="font-bold">RNA risk and initiation observations</h3>${rows}<p class="mt-2">Pattern hits and predicted MFE are computational observations, not measured expression or folding benefits. Missing observations are not PASS.</p></section>` : '';
}
function renderScientificReview(result) {
    let panel=document.getElementById('scientificReviewPanel');
    if (!panel) {
        panel=document.createElement('section');panel.id='scientificReviewPanel';
        elements.candidateComparisonContainer.after(panel);
    }
    panel.innerHTML=scientificReviewHtml({scientific_review:result.scientific_review,candidate_scientific_reviews:result.candidate_scientific_reviews || (result.candidates || []).map(c=>({id:c.id,label:c.label,scientific_review:c.scientific_review}))});
    panel.hidden=!result.scientific_review;
}

async function runOptimization() {
    // Flush the debounced input handler so an immediate click after paste/type
    // optimizes the current textarea value instead of stale state.
    const checkedDesignMode = document.querySelector('input[name="researchDesignMode"]:checked')?.value;
    if (checkedDesignMode) {
        state.researchDesignMode = checkedDesignMode;
    }
    handleSequenceChange({ target: elements.sequenceInput });

    if (!state.sequence || state.sequence.length < 3) {
        showToast('Please enter a valid DNA sequence', 'error');
        return;
    }

    if (state.host === 'human_stemcell') {
        showToast('Human Stem Cell (hMSC / iPSC) is currently In Development (research prototype). Online web optimization is locked pending wet-lab calibration. Please use factorforge.engines.stemcell_3dp locally.', 'info');
        return;
    }

    setLoading(true);
    trackEvent('optimization_run', {
        objective: state.objective,
        host: state.host,
        kozak: state.kozak,
        dinuc: state.dinuc,
        seq_len_bucket: seqLenBucket(state.sequence.length),
    });

    try {
        captureSopFromUi();
        // Prepare Request
        const payload = {
            sequence: state.sequence,
            host: state.host,
            mode: state.engineMode,
            save_db: state.saveDb,
            use_template: state.useTemplate,
            kozak: state.kozak,
            dinuc: state.dinuc,
            return_candidates: true,
            sop_profile: state.activeSop
        };
        if (state.engineMode !== 'profile') {
            payload.profile = 'balanced';
        } else if (['feasibility_best', 'dp_v2_1_1'].includes(state.objective) && state.host === 'nbenthamiana') {
            payload.objective = state.objective;
            payload.host_profile = state.host;
            // Host-aware default (v3.3.0) — previously hardcoded to
            // the legacy 55-65 band, which silently overrode the server's
            // resolve_host_gc_range() default for every feasibility_best run.
            payload.constraints = getGcRange(state.host);
        } else {
            payload.profile = ['feasibility_best', 'dp_v2_1_1'].includes(state.objective) ? 'balanced' : state.objective;
        }
        const seedValue = elements.optimizationSeed.value.trim();
        if (seedValue !== '') {
            const seed = Number(seedValue);
            if (!Number.isInteger(seed)) throw new Error('Seed must be an integer');
            payload.seed = seed;
        }
        const selectedTypeIisEnzymes = Array.from(elements.typeIisEnzymes)
            .filter(input => input.checked)
            .map(input => input.value);
        const customRestrictionSites = mergeRestrictionSitePresets(
            parseCustomRestrictionSites(elements.customRestrictionSites.value),
            selectedTypeIisEnzymes
        );
        state.customRestrictionSites = customRestrictionSites;
        state.selectedTypeIisEnzymes = selectedTypeIisEnzymes;
        if (customRestrictionSites.length > 0) {
            payload.custom_restriction_sites = customRestrictionSites;
        }
        state.structureViewerData = null;
        const rawFile = document.getElementById('structureWebFile')?.files?.[0] || null;
        const file = (rawFile && rawFile.size > 0 && rawFile.size <= 250000 && /\.(pdb|cif|mmcif)$/i.test(rawFile.name)) ? rawFile : null;
        const uploadedText = file ? await file.text() : (state.sampleStructureText || null);
        const uploadedFormat = file ? (/\.pdb$/i.test(file.name) ? 'pdb' : 'mmcif') : (state.sampleStructureFormat || 'pdb');
        const chainVal = document.getElementById('structureWebChain')?.value.trim() || state.sampleStructureChain || 'A';
        const modelVal = Number(document.getElementById('structureWebModel')?.value) || state.sampleStructureModel || 1;

        if (uploadedText) try {
            const resolveResp = await fetch(API_ENDPOINT, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    action: 'auto_resolve_structure',
                    sequence: state.sequence,
                    structure_context: uploadedText ? { content: uploadedText, format: uploadedFormat, chain: chainVal, model: modelVal } : null
                })
            });
            if (resolveResp.ok) {
                const resolveJson = await resolveResp.json();
                if (resolveJson?.resolved) {
                    state.structureViewerData = resolveJson.resolved;
                }
            }
        } catch (resolveErr) {
            console.warn('Auto structure resolution warning:', resolveErr);
        }

        if (state.researchDesignMode === 'paired_3dp') {
            if (apiCapabilities.direct_structure?.available !== true) throw new Error('3DP is unavailable on this backend. Update the API and UI together.');
            const resolved = state.structureViewerData;
            let rawPositions = document.getElementById('structureWebPositions').value.trim();
            if (!rawPositions && resolved?.positions?.length) {
                rawPositions = resolved.positions.join(',');
                document.getElementById('structureWebPositions').value = rawPositions;
            }
            if (!rawPositions) {
                throw new Error('Please specify pause positions (e.g. 105,106,107) or upload a PDB file to auto-detect them.');
            }
            let positions = [];
            const tokens = rawPositions.split(',').map(value => value.trim()).filter(Boolean);
            if (tokens.some(value => !/^\d+$/.test(value))) {
                throw new Error('Positions must be comma-separated residue numbers (e.g. 105,106,107)');
            }
            positions = tokens.map(Number);
            if (!positions.length || positions.some(p => !Number.isSafeInteger(p) || p < 1) || new Set(positions).size !== positions.length) {
                throw new Error('Positions must be unique residue numbers without duplicates (e.g. 105,106,107)');
            }
            const contentToUse = uploadedText || resolved?.content;
            const formatToUse = uploadedFormat || resolved?.format || 'pdb';
            if (!contentToUse) {
                throw new Error('Please upload an exact-match PDB/mmCIF structure file or click "Load Adalimumab (Humira) 3DP Sample".');
            }
            if (!document.getElementById('structureWebReviewed').checked) {
                document.getElementById('structureWebReviewed').checked = true;
            }
            if (resolved?.source_label) {
                document.getElementById('structureRecommendationStatus').textContent = `Structure source: ${resolved.source_label} · Pause positions (${positions.length}): ${positions.slice(0, 12).join(', ')}${positions.length > 12 ? '…' : ''}`;
            }
            payload.design_mode = 'paired_3dp';
            payload.mode = 'profile';
            payload.save_db = false;
            payload.kozak = false;
            payload.dinuc = false;
            payload.use_template = false;
            payload.constraints = getGcRange(state.host);
            delete payload.objective;
            delete payload.seed;
            payload.structure_context = {
                content: contentToUse,
                format: formatToUse,
                chain: chainVal,
                model: modelVal,
                positions,
                positions_reviewed: true,
                target_adaptiveness: Number(document.getElementById('structureWebTarget').value) || 0.40,
                penalty_weight: Number(document.getElementById('structureWebPenalty').value) || 2.0
            };
        }
        payload.acceptance_criteria = getAcceptanceCriteriaPayload();

        let response;
        try {
            response = await fetch(API_ENDPOINT, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(payload)
            });
        } catch (f) {
            console.warn('Network error or CORS issue', f);
            throw new Error('Network request failed');
        }

        let data;
        if (response.ok) {
            data = await response.json();
        } else {
            let message = `API request failed (${response.status})`;
            try {
                const errorBody = await response.json();
                if (errorBody && errorBody.error) message = errorBody.error;
            } catch (parseError) {
                console.warn('Unable to parse API error response', parseError);
            }
            throw new Error(message);
        }

        if (state.structureViewerData) {
            data.structure_viewer = {
                source_type: state.structureViewerData.source_type,
                source_label: state.structureViewerData.source_label,
                positions: state.structureViewerData.positions || [],
                ca_trace: state.structureViewerData.ca_trace || [],
                cartoon_pdb: state.structureViewerData.cartoon_pdb || null,
            };
        }
        data.input_sequence = state.sequence || '';
        state.results = data;
        addToHistory(state.sequence, data);
        renderResults();
        showToast('Optimization complete!', 'success');
        const primary = getPrimaryResult(data);
        if (primary?.metrics) {
            trackEvent('optimization_result', {
                objective: state.objective,
                host: getResultHostProfile(data),
                cai_bucket: caiBucket(primary.metrics.cai ?? 0),
                gc_bucket: gcBucket(primary.metrics.gc_percent ?? 0, getResultHostProfile(data)),
                success: true,
            });
        }

    } catch (error) {
        console.error('Optimization Failed:', error);
        elements.emptyState.classList.toggle('hidden', Boolean(state.results));
        elements.resultsContainer.classList.toggle('hidden', !state.results);

        if (ENABLE_MOCK && state.researchDesignMode !== 'paired_3dp') {
            showToast('Running mock optimization for demonstration...', 'info');
            await new Promise(r => setTimeout(r, 1500));
            state.results = getMockResult();
            renderResults();
            showToast('Showing simulated results', 'success');
        } else {
            if (ENABLE_CUSTOM_SOP) {
                document.getElementById('manualSopOverrides').open = true;
                document.getElementById('profileManagement').open = true;
            }
            showToast(`Optimization failed: ${error.message}`, 'error');
        }
    } finally {
        setLoading(false);
    }
}

// UI Updating Functions
function setLoading(loading) {
    state.isOptimizing = loading;
    elements.optimizeBtn.disabled = loading;
    if (loading) document.getElementById('settingsChangeNotice')?.classList.add('hidden');
    else handleSequenceChange({ target: elements.sequenceInput });

    if (typeof updateWorkflowStepper === 'function') {
        const currentView = elements.designWorkspace ? elements.designWorkspace.getAttribute('data-view-state') : 'configure';
        updateWorkflowStepper(loading ? 'review' : currentView, loading);
    }

    if (loading) {
        elements.btnText.classList.add('opacity-0');
        elements.loadingIndicator.classList.remove('hidden');
        elements.emptyState.classList.add('hidden');
        elements.resultsContainer.classList.add('hidden');
    } else {
        elements.btnText.classList.remove('opacity-0');
        elements.loadingIndicator.classList.add('hidden');
        elements.validationStatus.classList.remove('hidden');
    }
}

function updateStructureLinks() {
    const open3dBtn = document.getElementById('openLocal3dViewerBtn');
    if (open3dBtn && !open3dBtn.dataset.bound) {
        open3dBtn.addEventListener('click', () => {
            const container = document.getElementById('structure3dViewerContainer');
            if (container && !container.classList.contains('hidden')) {
                container.scrollIntoView({ behavior: 'smooth', block: 'center' });
            }
        });
        open3dBtn.dataset.bound = 'true';
    }
    const seq = (state.sequence || '').replace(/[^A-Za-z]/g, '').toUpperCase();
    if (!seq || seq.length < 10) {
        if (elements.alphafoldLink) elements.alphafoldLink.href = '#';
        if (elements.esmatlasFoldLink) elements.esmatlasFoldLink.href = '#';
        return;
    }
    const encoded = encodeURIComponent(seq);
    if (elements.alphafoldLink) {
        elements.alphafoldLink.href = `https://alphafold.ebi.ac.uk/search/sequence/${encoded}`;
    }
    if (elements.esmatlasFoldLink) {
        elements.esmatlasFoldLink.href = `https://esmatlas.com/resources?action=fold&sequence=${encoded}`;
    }
}

// Third-party structure-prediction linkout consent. AlphaFold DB and ESM Atlas
// have no data processing agreement with Eijex — this is a user-initiated
// transfer, not a subprocessor relationship, so the gate names the actual
// operator instead of implying an Eijex-controlled data flow.
const LINKOUT_CONSENT_COPY = {
    alphafold: {
        title: 'Open AlphaFold DB?',
        body: 'This will send your sequence to alphafold.ebi.ac.uk, operated by EMBL-EBI in partnership with Google DeepMind. Eijex does not control their logging, retention, or use of this data. Do not proceed if this sequence is confidential, proprietary, or controlled.'
    },
    esmatlas: {
        title: 'Open ESM Atlas?',
        body: 'This will send your sequence to esmatlas.com, operated by Meta Platforms, Inc. (Meta AI). Eijex does not control their logging, retention, or use of this data. Do not proceed if this sequence is confidential, proprietary, or controlled.'
    }
};

function openLinkoutConsent(serviceKey, targetHref) {
    const copy = LINKOUT_CONSENT_COPY[serviceKey];
    if (!copy || !targetHref) return;
    elements.linkoutConsentTitle.textContent = copy.title;
    elements.linkoutConsentBody.textContent = copy.body;
    elements.linkoutConsentContinue.href = targetHref;
    elements.linkoutConsentModal.classList.remove('hidden');
}

function closeLinkoutConsent() {
    elements.linkoutConsentModal.classList.add('hidden');
}

function getResultHostProfile(res) {
    return res?.host_profile || res?.validation?.host_profile || state.host || 'nbenthamiana';
}

function formatHostProfile(hostProfile) {
    const normalized = String(hostProfile || 'nbenthamiana').toLowerCase();
    const label = HOST_LABELS[normalized] || hostProfile;
    return `${hostProfile} (${label})`;
}

function normalizedDna(value) {
    return String(value || '').toUpperCase().replace(/[^ACGT]/g, '');
}

function codonsWithoutTerminalStop(sequence) {
    const dna = normalizedDna(sequence);
    if (!dna || dna.length % 3 !== 0) return { valid: false, codons: [], hasStopCodon: false };
    const codons = dna.match(/.{3}/g) || [];
    const hasStopCodon = CODON_TABLE[codons.at(-1)] === '*';
    return { valid: true, codons: hasStopCodon ? codons.slice(0, -1) : codons, hasStopCodon };
}

function responseCodonFrequencies(res, hostId) {
    const candidates = [
        res?.host_codon_frequencies,
        res?.codon_frequencies,
        res?.host_metadata?.[hostId]?.codon_frequencies,
        hostMetadata?.[hostId]?.codon_frequencies,
    ];
    const source = candidates.find(value => value && typeof value === 'object');
    if (!source) return {};
    return Object.fromEntries(Object.entries(source).map(([codon, raw]) => {
        const value = Number(typeof raw === 'object' ? raw.frequency : raw);
        return [codon.toUpperCase(), Number.isFinite(value) ? (value > 1 ? value / 100 : value) : null];
    }));
}

function localGcAtCodon(sequence, codonIndex, windowSize = 50) {
    const dna = normalizedDna(sequence);
    if (!dna) return null;
    const center = codonIndex * 3 + 1;
    const start = Math.max(0, Math.min(dna.length - windowSize, center - Math.floor(windowSize / 2)));
    const window = dna.slice(start, start + windowSize);
    return window ? ((window.match(/[GC]/g) || []).length / window.length) * 100 : null;
}

function buildCodonAccounting(originalCds, optimizedCds, res, gcTarget) {
    const original = codonsWithoutTerminalStop(originalCds);
    const optimized = codonsWithoutTerminalStop(optimizedCds);
    const comparable = original.valid && optimized.valid && original.codons.length > 0;
    if (!comparable) return { available: false, reason: originalCds ? 'Valid CDS reference required' : 'CDS reference required', codons: [] };
    const frequencies = responseCodonFrequencies(res, getResultHostProfile(res));
    const length = Math.max(original.codons.length, optimized.codons.length);
    const codons = [];
    let synonymousChanges = 0, unchanged = 0, substitutions = 0, insertions = 0, deletions = 0;
    let higherPref = 0, lowerPref = 0, neutralPref = 0, unknownPref = 0;
    for (let index = 0; index < length; index += 1) {
        const origCodon = original.codons[index] || null;
        const optCodon = optimized.codons[index] || null;
        const origAa = origCodon ? CODON_TABLE[origCodon] : null;
        const optAa = optCodon ? CODON_TABLE[optCodon] : null;
        if (!origCodon) insertions += 1;
        else if (!optCodon) deletions += 1;
        else if (origAa !== optAa) substitutions += 1;
        else if (origCodon === optCodon) unchanged += 1;
        else synonymousChanges += 1;
        const origFreq = origCodon && Object.hasOwn(frequencies, origCodon) ? frequencies[origCodon] : null;
        const optFreq = optCodon && Object.hasOwn(frequencies, optCodon) ? frequencies[optCodon] : null;
        if (origCodon && optCodon && origCodon !== optCodon && origAa === optAa) {
            if (origFreq == null || optFreq == null) unknownPref += 1;
            else if (optFreq > origFreq) higherPref += 1;
            else if (optFreq < origFreq) lowerPref += 1;
            else neutralPref += 1;
        }
        codons.push({
            index, aa: optAa || origAa || '?', origAa, optAa, origCodon, optCodon,
            isChanged: origCodon !== optCodon, origFreq, optFreq,
            localGc: localGcAtCodon(optimizedCds, index),
        });
    }
    const totalCodons = original.codons.length;
    const comparedResidues = Math.min(original.codons.length, optimized.codons.length);
    const aaIdentity = totalCodons ? (Math.max(0, comparedResidues - substitutions) / Math.max(original.codons.length, optimized.codons.length)) * 100 : 0;
    return {
        available: true, totalCodons, synonymousChanges, unchanged, substitutions, insertions, deletions, aaIdentity,
        higherPref, lowerPref, neutralPref, unknownPref, hasStopCodon: original.hasStopCodon,
        aaPreserved: substitutions === 0 && insertions === 0 && deletions === 0,
        gcTarget, codons,
    };
}

function percentOf(value, total) {
    return total ? `${((value / total) * 100).toFixed(1)}%` : '0.0%';
}

function renderCodonAccounting(accounting) {
    if (!accounting.available) {
        elements.statTotalCodons.textContent = 'N/A'; elements.statShifts.textContent = 'N/A';
        elements.statUnchanged.textContent = 'N/A'; elements.statSubstitutions.textContent = 'N/A';
        elements.statShiftBreakdown.textContent = accounting.reason;
        elements.statAaIdentity.textContent = 'Translation comparison unavailable';
        elements.trackAvailability.textContent = accounting.reason;
        return;
    }
    elements.statTotalCodons.textContent = String(accounting.totalCodons);
    elements.statShifts.textContent = `${accounting.synonymousChanges} / ${percentOf(accounting.synonymousChanges, accounting.totalCodons)}`;
    elements.statUnchanged.textContent = `${accounting.unchanged} / ${percentOf(accounting.unchanged, accounting.totalCodons)}`;
    elements.statSubstitutions.textContent = `${accounting.substitutions} substitutions`;
    const known = `Higher ${accounting.higherPref} · Lower ${accounting.lowerPref} · Neutral ${accounting.neutralPref}`;
    elements.statShiftBreakdown.textContent = accounting.unknownPref ? `${known} · Frequency N/A ${accounting.unknownPref}` : known;
    elements.statAaIdentity.textContent = `${accounting.aaIdentity.toFixed(1)}% identity · ${accounting.insertions} insertions · ${accounting.deletions} deletions`;
    elements.trackAvailability.textContent = `${accounting.totalCodons} sense codons · terminal stop ${accounting.hasStopCodon ? 'excluded' : 'not present'}`;
}

function frequencyColor(frequency) {
    if (frequency == null) return '#334155';
    const percent = frequency * 100;
    if (percent < 20) return '#ef4444'; if (percent < 35) return '#f59e0b';
    if (percent < 50) return '#64748b'; if (percent <= 70) return '#22c55e'; return '#10b981';
}

function drawCodonTrack(canvas, codons, frequencyKey, hoverIndex = null) {
    if (!canvas) return;
    const rect = canvas.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    const width = Math.max(1, rect.width); const height = 56;
    canvas.width = Math.round(width * dpr); canvas.height = Math.round(height * dpr);
    const ctx = canvas.getContext('2d'); ctx.scale(dpr, dpr); ctx.clearRect(0, 0, width, height);
    ctx.fillStyle = '#020617'; ctx.fillRect(0, 0, width, height);
    if (!codons.length) { ctx.fillStyle = '#64748b'; ctx.font = '11px sans-serif'; ctx.fillText('No comparable CDS track', 12, 32); return; }
    const segmentWidth = width / codons.length;
    codons.forEach((codon, index) => {
        ctx.fillStyle = frequencyColor(codon[frequencyKey]);
        ctx.fillRect(index * segmentWidth, 13, Math.max(1, segmentWidth + .25), 30);
    });
    if (hoverIndex != null) {
        const x = Math.min(width - .5, (hoverIndex + .5) * segmentWidth);
        ctx.strokeStyle = '#f8fafc'; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(x, 2); ctx.lineTo(x, 54); ctx.stroke();
    }
}

function formatFrequency(value) { return value == null ? 'N/A' : `${(value * 100).toFixed(1)}%`; }

function updateCodonInspector(accounting, index) {
    const item = accounting.codons[index]; if (!item) return;
    const name = AA_NAMES[item.aa] || item.aa;
    elements.inspectorPosition.textContent = `Codon ${index + 1} / ${name} (${item.aa}) · nt ${index * 3 + 1}–${index * 3 + 3}`;
    elements.inspectorCodons.textContent = `${item.origCodon || '—'} (${formatFrequency(item.origFreq)}) → ${item.optCodon || '—'} (${formatFrequency(item.optFreq)})`;
    const effects = [];
    if (!item.isChanged) effects.push('Codon unchanged');
    else if (item.origAa === item.optAa) effects.push('Synonymous codon shift observed');
    else effects.push(`Translated residue changed: ${item.origAa || '—'} → ${item.optAa || '—'}`);
    if (item.decision_rationale) {
        effects.push(`Trace: ${item.decision_rationale}`);
    } else {
        // Strict Rule: No UI inference
    }
    elements.inspectorRationale.textContent = effects.join(' · ');
    elements.inspectorConstraints.textContent = accounting.aaPreserved ? 'PASS · amino-acid sequence preserved' : 'FAIL · amino-acid sequence differs';
}

function renderCodonTracks(accounting) {
    const codons = accounting.available ? accounting.codons : [];
    if (codons.length > 0) {
        elements.canvasOriginalTrack?.classList.add('has-track');
        elements.canvasOptimizedTrack?.classList.add('has-track');
    } else {
        elements.canvasOriginalTrack?.classList.remove('has-track');
        elements.canvasOptimizedTrack?.classList.remove('has-track');
    }
    const draw = index => { drawCodonTrack(elements.canvasOriginalTrack, codons, 'origFreq', index); drawCodonTrack(elements.canvasOptimizedTrack, codons, 'optFreq', index); };
    draw(null);
    if (!elements.interactiveTrackContainer.dataset.bound) {
        elements.interactiveTrackContainer.addEventListener('mousemove', event => {
            const current = state.reportData?.accounting;
            if (!current?.available || !current.codons.length) return;
            const rect = elements.canvasOptimizedTrack.getBoundingClientRect();
            const index = Math.max(0, Math.min(current.codons.length - 1, Math.floor(((event.clientX - rect.left) / rect.width) * current.codons.length)));
            drawCodonTrack(elements.canvasOriginalTrack, current.codons, 'origFreq', index);
            drawCodonTrack(elements.canvasOptimizedTrack, current.codons, 'optFreq', index);
            updateCodonInspector(current, index);
        });
        elements.interactiveTrackContainer.addEventListener('mouseleave', () => renderCodonTracks(state.reportData?.accounting || { available:false, codons:[] }));
        window.addEventListener('resize', () => renderCodonTracks(state.reportData?.accounting || { available:false, codons:[] }));
        elements.interactiveTrackContainer.dataset.bound = 'true';
    }
    if (accounting.available) updateCodonInspector(accounting, 0);
}

const FALLBACK_NB_CODON_FREQ = Object.freeze({
    TTT:0.58,TTC:0.42,TTA:0.14,TTG:0.24,CTT:0.26,CTC:0.14,CTA:0.09,CTG:0.13,
    ATT:0.51,ATC:0.29,ATA:0.20,ATG:1.00,GTT:0.41,GTC:0.18,GTA:0.15,GTG:0.26,
    TCT:0.26,TCC:0.14,TCA:0.22,TCG:0.07,AGT:0.17,AGC:0.14,CCT:0.38,CCC:0.13,
    CCA:0.38,CCG:0.11,ACT:0.40,ACC:0.20,ACA:0.31,ACG:0.09,GCT:0.44,GCC:0.17,
    GCA:0.31,GCG:0.08,TAT:0.57,TAC:0.43,CAT:0.61,CAC:0.39,CAA:0.56,CAG:0.44,
    AAT:0.59,AAC:0.41,AAA:0.47,AAG:0.53,GAT:0.67,GAC:0.33,GAA:0.55,GAG:0.45,
    TGT:0.59,TGC:0.41,TGG:1.00,CGT:0.16,CGC:0.07,CGA:0.10,CGG:0.07,AGA:0.35,
    AGG:0.25,GGT:0.34,GGC:0.15,GGA:0.36,GGG:0.15,TAA:0.45,TAG:0.20,TGA:0.35
});

function buildFallbackCaTrace(codons) {
    const n = codons.length;
    const trace = [];
    for (let i = 0; i < n; i++) {
        const codon = codons[i] || 'GCT';
        const aa = CODON_TABLE[codon] || 'A';
        const seg = Math.floor(i / 28);
        const idxInSeg = i % 28;
        const xBase = seg * 12.0;
        const direction = seg % 2 === 0 ? 1.0 : -1.0;
        let x, y, z;
        if (idxInSeg < 20) {
            const theta = idxInSeg * (100.0 * Math.PI / 180.0);
            x = xBase + 2.3 * Math.cos(theta);
            y = 2.3 * Math.sin(theta);
            z = direction * (idxInSeg * 1.5 - 15.0);
        } else {
            const loopT = (idxInSeg - 19) / 9.0;
            const nextXBase = (seg + 1) * 12.0;
            x = xBase * (1.0 - loopT) + nextXBase * loopT;
            y = 4.0 * Math.sin(loopT * Math.PI) * direction;
            z = direction * (15.0 + 4.5 * Math.sin(loopT * Math.PI));
        }
        trace.push({
            resSeq: i + 1,
            resName: (AA_NAMES[aa] || 'ALA').toUpperCase(),
            aa,
            x: Number(x.toFixed(2)),
            y: Number(y.toFixed(2)),
            z: Number(z.toFixed(2)),
            bFactor: 85.0
        });
    }
    return trace;
}

function synthesizeBaselineReferenceDna(optSeq, res) {
    const parsed = codonsWithoutTerminalStop(optSeq || '');
    if (!parsed.valid || !parsed.codons.length) return '';
    const freqsRaw = responseCodonFrequencies(res, getResultHostProfile(res));
    const freqs = Object.keys(freqsRaw).length ? freqsRaw : FALLBACK_NB_CODON_FREQ;
    const is3dp = res?.design_mode === 'paired_3dp' || res?.profile === 'structure_3dp';
    const byAa = {};
    for (const [codon, aa] of Object.entries(CODON_TABLE)) {
        if (aa === '*') continue;
        if (!byAa[aa]) byAa[aa] = [];
        byAa[aa].push({ codon, f: freqs[codon] ?? FALLBACK_NB_CODON_FREQ[codon] ?? 0.2 });
    }
    for (const aa of Object.keys(byAa)) {
        byAa[aa].sort((a, b) => b.f - a.f);
    }
    const refCodons = parsed.codons.map(c => {
        const aa = CODON_TABLE[c];
        const list = byAa[aa] || [];
        if (!list.length) return c;
        if (is3dp) {
            // 1DP baseline uses the top-frequency synonymous codon
            return list[0].codon;
        }
        // For 1DP protein-only input, compare against canonical uniform median codon
        return list.length > 1 ? list[Math.min(1, list.length - 1)].codon : list[0].codon;
    });
    const stopCodon = (optSeq && optSeq.length >= parsed.codons.length * 3 + 3)
        ? optSeq.slice(parsed.codons.length * 3, parsed.codons.length * 3 + 3)
        : 'TAA';
    return refCodons.join('') + stopCodon;
}

function buildStructureViewerPayload(res, primary) {
    const is3dp = res?.design_mode === 'paired_3dp' || res?.profile === 'structure_3dp';
    const freqsRaw = responseCodonFrequencies(res, getResultHostProfile(res));
    const freqs = Object.keys(freqsRaw).length ? freqsRaw : FALLBACK_NB_CODON_FREQ;
    const candidateList = Array.isArray(res?.candidates)
        ? res.candidates
        : (Array.isArray(res?.report_candidates) ? res.report_candidates : []);
    const cand1dp = candidateList.find(c => c.id === 'baseline_1dp') || null;
    const cand3dp = candidateList.find(c => c.id === 'structure_3dp') || null;
    const rawInputSeq = res?.input_sequence || state.sequence || '';
    const seqOrigRaw = !isProteinInputResult(res) ? rawInputSeq : '';
    const optPrimarySeq = primary?.optimized_sequence || '';
    const seq1dp = is3dp
        ? (cand1dp?.dna_sequence || synthesizeBaselineReferenceDna(optPrimarySeq, res) || optPrimarySeq)
        : optPrimarySeq;
    const seq3dp = is3dp ? (cand3dp?.dna_sequence || optPrimarySeq) : null;
    const parsedOrig = codonsWithoutTerminalStop(seqOrigRaw);
    const parsed1dp = codonsWithoutTerminalStop(seq1dp);
    const parsed3dp = seq3dp ? codonsWithoutTerminalStop(seq3dp) : { valid: false, codons: [] };
    const codonsOrig = parsedOrig.valid ? parsedOrig.codons : [];
    const codons1dp = parsed1dp.codons || [];
    const codons3dp = parsed3dp.codons || [];
    const total = Math.max(codons1dp.length, codons3dp.length);
    const storedViewer = res?.structure_viewer || null;
    const expectedProtein = codons1dp.map(c => CODON_TABLE[c] || 'X').join('');
    const matchesProtein = storedViewer?.ca_trace?.length === total
        && storedViewer.ca_trace.map(p => p.aa).join('') === expectedProtein;
    if (!matchesProtein) return {mode: is3dp ? '3dp' : '1dp', source_label: 'NOT_AVAILABLE', cartoonPdb: null, residues: []};
    const caTrace = storedViewer.ca_trace;
    const pauseSet = new Set(storedViewer?.positions || []);
    const residues = [];
    for (let i = 0; i < total; i++) {
        const c1 = codons1dp[i] || codons3dp[i] || 'GCT';
        const c3 = is3dp ? (codons3dp[i] || c1) : null;
        const cOrig = codonsOrig[i] || c1;
        const fOrig = freqs[cOrig] ?? FALLBACK_NB_CODON_FREQ[cOrig] ?? 0.35;
        const f1 = freqs[c1] ?? FALLBACK_NB_CODON_FREQ[c1] ?? 0.35;
        const f3 = c3 ? (freqs[c3] ?? FALLBACK_NB_CODON_FREQ[c3] ?? 0.35) : null;
        const pt = caTrace[i] || { resSeq: i + 1, aa: CODON_TABLE[c1] || 'A', x: i * 1.5, y: 0, z: 0, plddt: 85.0 };
        const plddtVal = Number(pt.plddt ?? pt.bFactor ?? 85.0);
        residues.push({
            pos: i + 1,
            aa: CODON_TABLE[c3 || c1] || pt.aa || 'A',
            x: Number(Array.isArray(pt.ca) ? pt.ca[0] : (pt.x || 0)),
            y: Number(Array.isArray(pt.ca) ? pt.ca[1] : (pt.y || 0)),
            z: Number(Array.isArray(pt.ca) ? pt.ca[2] : (pt.z || 0)),
            plddt: Number(plddtVal.toFixed(1)),
            codonOrig: cOrig,
            freqOrig: Number(fOrig.toFixed(4)),
            hasOrigDna: Boolean(codonsOrig.length),
            codon1dp: c1,
            freq1dp: Number(f1.toFixed(4)),
            codon3dp: c3,
            freq3dp: f3 != null ? Number(f3.toFixed(4)) : null,
            isDiffOrig1dp: Boolean(codonsOrig.length ? cOrig !== c1 : f1 >= 0.35),
            isDiffOrig3dp: Boolean(is3dp && (codonsOrig.length ? cOrig !== c3 : c1 !== c3)),
            isDiff: Boolean(is3dp && c1 !== c3),
            isPause: pauseSet.has(i + 1)
        });
    }
    const defaultPauseFocus = residues.find(r => r.isPause && r.isDiff) || residues.find(r => r.isPause) || residues[Math.floor(residues.length / 2)] || residues[0];
    return {
        mode: is3dp ? '3dp' : '1dp',
        source_label: matchesProtein ? (storedViewer?.source_label || 'Uploaded structure (provider unverified)') : 'NOT_AVAILABLE',
        cartoonPdb: matchesProtein ? (storedViewer?.cartoon_pdb || null) : null,
        cai1dp: Number(cand1dp?.cai ?? primary?.metrics?.cai ?? 0),
        gc1dp: Number(cand1dp?.gc_percent ?? primary?.metrics?.gc_percent ?? 0),
        cai3dp: is3dp ? Number(cand3dp?.cai ?? primary?.metrics?.cai ?? 0) : null,
        gc3dp: is3dp ? Number(cand3dp?.gc_percent ?? primary?.metrics?.gc_percent ?? 0) : null,
        delta3dGain: is3dp && res?.structure_comparison ? Number(res.structure_comparison.delta_3d_gain || 0) : null,
        totalEdits: is3dp && res?.structure_comparison ? Number(res.structure_comparison.total_edits || 0) : residues.filter(r => r.isDiff).length,
        pauseCount: pauseSet.size,
        defaultSelectedPos: defaultPauseFocus ? defaultPauseFocus.pos : 1,
        residues
    };
}

function ff3dCatmullRom(p0, p1, p2, p3, t) {
    const t2 = t * t, t3 = t2 * t;
    return [
        0.5 * ((2 * p1[0]) + (-p0[0] + p2[0]) * t + (2 * p0[0] - 5 * p1[0] + 4 * p2[0] - p3[0]) * t2 + (-p0[0] + 3 * p1[0] - 3 * p2[0] + p3[0]) * t3),
        0.5 * ((2 * p1[1]) + (-p0[1] + p2[1]) * t + (2 * p0[1] - 5 * p1[1] + 4 * p2[1] - p3[1]) * t2 + (-p0[1] + 3 * p1[1] - 3 * p2[1] + p3[1]) * t3),
        0.5 * ((2 * p1[2]) + (-p0[2] + p2[2]) * t + (2 * p0[2] - 5 * p1[2] + 4 * p2[2] - p3[2]) * t2 + (-p0[2] + 3 * p1[2] - 3 * p2[2] + p3[2]) * t3)
    ];
}

function ff3dShadeColor(hex, factor) {
    const h = String(hex || '#38bdf8').replace('#', '');
    const r = Math.min(255, Math.max(0, Math.round(parseInt(h.substring(0, 2), 16) * factor)));
    const g = Math.min(255, Math.max(0, Math.round(parseInt(h.substring(2, 4), 16) * factor)));
    const b = Math.min(255, Math.max(0, Math.round(parseInt(h.substring(4, 6), 16) * factor)));
    return `rgb(${r},${g},${b})`;
}

// Standard AlphaFold Server pLDDT color scale
function ff3dColorLeftInput(r) {
    if (r.plddt >= 90) return '#0053d6';
    if (r.plddt >= 70) return '#65cbf3';
    if (r.plddt >= 50) return '#ffdb13';
    return '#ff7d45';
}

function ff3dColorRightPost(r, modeKey) {
    if (modeKey === '3dp') {
        if (r.isPause && r.isDiff) return '#ef4444';
        if (r.isPause) return '#f59e0b';
        if (r.isDiff) return '#10b981';
        return '#65cbf3';
    }
    if (r.isDiffOrig1dp) {
        return r.freq1dp >= 0.42 ? '#10b981' : '#0053d6';
    }
    return '#65cbf3';
}

function attachInteractive3dCanvas(canvas, residues, getModeKey, getCam, setCam, getSelectedPos, onSelectPos, onHover) {
    if (!canvas || !residues?.length) return () => {};
    let dragging = false, moved = false, lx = 0, ly = 0;
    const cx = residues.reduce((s, r) => s + r.x, 0) / residues.length;
    const cy = residues.reduce((s, r) => s + r.y, 0) / residues.length;
    const cz = residues.reduce((s, r) => s + r.z, 0) / residues.length;
    const pts = residues.map(r => ({ ...r, cx: r.x - cx, cy: r.y - cy, cz: r.z - cz }));
    let projectedResidues = [];

    function rotateAndProject(pt3, w, h, cam) {
        const cosY = Math.cos(cam.ry), sinY = Math.sin(cam.ry);
        const cosX = Math.cos(cam.rx), sinX = Math.sin(cam.rx);
        const x1 = pt3[0] * cosY + pt3[2] * sinY;
        const z1 = -pt3[0] * sinY + pt3[2] * cosY;
        const y2 = pt3[1] * cosX - z1 * sinX;
        const z2 = pt3[1] * sinX + z1 * cosX;
        const fov = 560 * cam.zoom;
        const dist = 160;
        const scale = fov / (dist + z2);
        return {
            sx: w * 0.5 + x1 * scale + (cam.panX || 0),
            sy: h * 0.5 - y2 * scale + (cam.panY || 0),
            z: z2,
            scale
        };
    }

    function draw() {
        const rect = canvas.getBoundingClientRect();
        if (!rect.width || !rect.height) return;
        const dpr = window.devicePixelRatio || 1;
        canvas.width = Math.round(rect.width * dpr);
        canvas.height = Math.round(rect.height * dpr);
        const ctx = canvas.getContext('2d');
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        const w = rect.width, h = rect.height;

        ctx.fillStyle = '#f8fafc';
        ctx.fillRect(0, 0, w, h);

        const cam = getCam();
        const modeKey = typeof getModeKey === 'function' ? getModeKey() : getModeKey;
        const selectedPos = typeof getSelectedPos === 'function' ? getSelectedPos() : null;
        const steps = 4;
        const segments = [];
        for (let i = 0; i < pts.length - 1; i++) {
            const p0 = pts[Math.max(0, i - 1)];
            const p1 = pts[i];
            const p2 = pts[i + 1];
            const p3 = pts[Math.min(pts.length - 1, i + 2)];
            const c0 = [p0.cx, p0.cy, p0.cz], c1 = [p1.cx, p1.cy, p1.cz], c2 = [p2.cx, p2.cy, p2.cz], c3 = [p3.cx, p3.cy, p3.cz];
            const baseCol = modeKey === 'orig' ? ff3dColorLeftInput(p1) : ff3dColorRightPost(p1, modeKey);
            const isHi = modeKey === 'orig'
                ? true
                : (modeKey === '3dp' ? (p1.isPause || p1.isDiff) : p1.isDiffOrig1dp);
            for (let s = 0; s < steps; s++) {
                const tA = s / steps, tB = (s + 1) / steps;
                const ptA = rotateAndProject(ff3dCatmullRom(c0, c1, c2, c3, tA), w, h, cam);
                const ptB = rotateAndProject(ff3dCatmullRom(c0, c1, c2, c3, tB), w, h, cam);
                segments.push({
                    type: 'seg',
                    z: (ptA.z + ptB.z) * 0.5,
                    ptA,
                    ptB,
                    color: baseCol,
                    width: (isHi ? 2.5 : 1.9) * ((ptA.scale + ptB.scale) * 0.5),
                    alpha: 0.96
                });
            }
        }

        projectedResidues = pts.map(r => {
            const p = rotateAndProject([r.cx, r.cy, r.cz], w, h, cam);
            return { ...p, r };
        });

        for (const pr of projectedResidues) {
            const r = pr.r;
            const isSel = r.pos === selectedPos;
            const showSphere = isSel || (modeKey === '3dp' && r.isPause && r.isDiff);
            if (showSphere) {
                segments.push({
                    type: 'node',
                    z: pr.z - (isSel ? 5 : 0),
                    pr,
                    isSel
                });
            }
        }

        segments.sort((a, b) => b.z - a.z);

        for (const item of segments) {
            const depthFactor = Math.max(0.55, Math.min(1.18, 1.0 - item.z / 110));
            if (item.type === 'seg') {
                ctx.save();
                ctx.globalAlpha = item.alpha;
                ctx.beginPath();
                ctx.moveTo(item.ptA.sx, item.ptA.sy);
                ctx.lineTo(item.ptB.sx, item.ptB.sy);
                ctx.strokeStyle = ff3dShadeColor(item.color, depthFactor);
                ctx.lineWidth = Math.max(2.0, item.width);
                ctx.lineCap = 'round';
                ctx.stroke();
                ctx.restore();
            } else if (item.type === 'node') {
                const { pr, isSel } = item;
                const radius = Math.max(3.0, (isSel ? 2.8 : 1.7) * pr.scale);
                ctx.save();
                ctx.beginPath();
                ctx.arc(pr.sx, pr.sy, radius, 0, Math.PI * 2);
                ctx.fillStyle = isSel ? '#0f172a' : '#ef4444';
                ctx.fill();
                ctx.lineWidth = 2.0;
                ctx.strokeStyle = '#ffffff';
                ctx.stroke();
                ctx.restore();
            }
        }
    }

    canvas.addEventListener('mousedown', e => {
        dragging = true;
        moved = false;
        lx = e.clientX;
        ly = e.clientY;
    });
    window.addEventListener('mouseup', () => { dragging = false; });
    window.addEventListener('mousemove', e => {
        if (!dragging) return;
        const dx = e.clientX - lx, dy = e.clientY - ly;
        if (Math.abs(dx) > 2 || Math.abs(dy) > 2) moved = true;
        const cam = getCam();
        if (e.shiftKey) {
            setCam({ ...cam, panX: (cam.panX || 0) + dx, panY: (cam.panY || 0) + dy });
        } else {
            setCam({
                ...cam,
                rx: Math.max(-1.45, Math.min(1.45, cam.rx + dy * 0.008)),
                ry: cam.ry + dx * 0.008
            });
        }
        lx = e.clientX;
        ly = e.clientY;
    });
    canvas.addEventListener('click', e => {
        if (moved) return;
        const rect = canvas.getBoundingClientRect();
        const mx = e.clientX - rect.left, my = e.clientY - rect.top;
        let best = null, minD = 18;
        for (const p of projectedResidues) {
            const d = Math.hypot(p.sx - mx, p.sy - my);
            if (d < minD) { minD = d; best = p; }
        }
        if (best && onSelectPos) onSelectPos(best.r.pos);
    });
    canvas.addEventListener('mousemove', e => {
        if (dragging) return;
        const rect = canvas.getBoundingClientRect();
        const mx = e.clientX - rect.left, my = e.clientY - rect.top;
        let best = null, minD = 16;
        for (const p of projectedResidues) {
            const d = Math.hypot(p.sx - mx, p.sy - my);
            if (d < minD) { minD = d; best = p; }
        }
        if (onHover) onHover(best ? best.r : null);
    });
    canvas.addEventListener('mouseleave', () => { if (onHover) onHover(null); });
    canvas.addEventListener('wheel', e => {
        e.preventDefault();
        const cam = getCam();
        const nextZoom = Math.max(0.45, Math.min(3.5, cam.zoom * (e.deltaY > 0 ? 0.92 : 1.08)));
        setCam({ ...cam, zoom: nextZoom });
    }, { passive: false });
    return draw;
}

function renderStructure3dViewer(res, primary) {
    const container = document.getElementById('structure3dViewerContainer');
    if (!container) return;
    const payload = buildStructureViewerPayload(res, primary);
    if (!payload.residues.length || !payload.cartoonPdb) {
        container.classList.add('hidden');
        return;
    }
    container.classList.remove('hidden');

    const is3dp = payload.mode === '3dp';
    let leftMode = 'orig';
    let rightMode = is3dp ? '3dp' : '1dp';
    let selectedPos = payload.defaultSelectedPos || 1;
    let syncCamera = false;
    let spin = false;
    let viewerBg = '#ffffff'; // Default AlphaFold Server white background

    const residueByPos = new Map(payload.residues.map(r => [r.pos, r]));

    const titleText = is3dp
        ? 'Original Structure Coordinates vs 3DP Structure-Aware Codon Mapping'
        : 'Original Structure Coordinates vs 1DP Baseline Codon Mapping';
    const subText = is3dp
        ? `Structure source: ${escapeHtml(payload.source_label)} · Comparing Left (pLDDT confidence) vs Right (3DP translational pause targets, ${payload.pauseCount} residues) on the identical 3D backbone.`
        : `Structure source: ${escapeHtml(payload.source_label)} · Comparing Left (pLDDT confidence) vs Right (1DP preferred-codon substitutions) on the identical 3D backbone.`;

    const modeButtonsHtml = is3dp ? `
        <div class="flex flex-wrap items-center gap-1.5">
            <button type="button" data-3d-pair="orig_3dp" class="px-3 py-1.5 rounded-lg text-[11px] font-extrabold bg-emerald-600 text-white shadow-sm">Original | 3DP</button>
            <button type="button" data-3d-pair="orig_1dp" class="px-3 py-1.5 rounded-lg text-[11px] font-bold bg-slate-800 text-slate-300 hover:bg-slate-700">Original | 1DP</button>
            <button type="button" data-3d-pair="1dp_3dp" class="px-3 py-1.5 rounded-lg text-[11px] font-bold bg-slate-800 text-slate-300 hover:bg-slate-700">1DP | 3DP</button>
        </div>
    ` : `<span class="px-2.5 py-1 rounded-lg bg-sky-900/60 border border-sky-500/40 text-sky-200 text-[10px] font-extrabold">Original | 1DP Side-by-Side</span>`;

    container.innerHTML = `
        <div class="codon-track-head">
            <div>
                <p class="text-[10px] font-extrabold uppercase tracking-widest text-emerald-400">${is3dp ? 'Uploaded Coordinate Ribbon · Original vs 3DP' : 'Uploaded Coordinate Ribbon · Original vs 1DP'}</p>
                <h3 class="text-base font-black text-white mt-0.5">${titleText}</h3>
                <p class="text-[11px] text-slate-300 mt-1">${subText}</p>
            </div>
            <div class="flex flex-wrap items-center gap-2">
                ${modeButtonsHtml}
                <button type="button" id="ff3dBgBtn" class="px-2.5 py-1.5 rounded-lg bg-slate-800 hover:bg-slate-700 text-slate-200 text-[10px] font-bold border border-slate-700">Background: White (AlphaFold)</button>
                <button type="button" id="ff3dSyncBtn" class="px-2.5 py-1.5 rounded-lg bg-slate-800 hover:bg-slate-700 text-slate-200 text-[10px] font-bold border border-slate-700">Camera Sync: Off (Independent)</button>
                <button type="button" id="ff3dSpinBtn" class="px-2.5 py-1.5 rounded-lg bg-slate-800 hover:bg-slate-700 text-slate-200 text-[10px] font-bold border border-slate-700">Auto Spin: Off</button>
                <button type="button" id="ff3dResetBtn" class="px-2.5 py-1.5 rounded-lg bg-slate-800 hover:bg-slate-700 text-slate-200 text-[10px] font-bold border border-slate-700">Reset View</button>
            </div>
        </div>
        <div id="ff3dGridWrap" class="grid grid-cols-1 lg:grid-cols-2 gap-3 mt-3">
            <div id="ff3dPaneLeft" class="relative rounded-xl overflow-hidden border border-slate-700 bg-white" style="height:430px;">
                <div class="absolute top-3 left-3 right-3 z-20 pointer-events-none flex items-start justify-between gap-2 bg-slate-900/90 border border-slate-700/80 rounded-lg px-3 py-2 backdrop-blur-sm shadow-sm">
                    <div>
                        <p id="ff3dLeftTitle" class="text-xs font-extrabold text-sky-400">LEFT · Original Structure (Input Structure)</p>
                        <p id="ff3dLeftSub" class="text-[10px] text-slate-300">Input coordinates · 3D cartoon ribbon (pLDDT confidence colors)</p>
                    </div>
                    <span id="ff3dLeftBadge" class="px-2 py-0.5 rounded-full text-[10px] font-extrabold bg-sky-950 text-sky-300 border border-sky-500/40 shrink-0">Input coordinates</span>
                </div>
                <div id="ff3dViewerLeft" class="w-full h-full relative"><canvas id="ff3dCanvasLeft" class="w-full h-full block cursor-grab"></canvas></div>
                <div id="ff3dLeftLegend" class="absolute bottom-2.5 left-2.5 right-2.5 z-20 pointer-events-none flex flex-wrap items-center gap-3 bg-slate-900/90 border border-slate-700 rounded-lg px-2.5 py-1.5 text-[10px] text-slate-200 shadow-sm">
                    <span class="inline-flex items-center gap-1"><span class="w-2.5 h-2.5 rounded-sm" style="background:#0053d6"></span>Very High (pLDDT ≥ 90)</span>
                    <span class="inline-flex items-center gap-1"><span class="w-2.5 h-2.5 rounded-sm" style="background:#65cbf3"></span>Confident (90 &gt; pLDDT ≥ 70)</span>
                    <span class="inline-flex items-center gap-1"><span class="w-2.5 h-2.5 rounded-sm" style="background:#ffdb13"></span>Low (70 &gt; pLDDT ≥ 50)</span>
                    <span class="inline-flex items-center gap-1"><span class="w-2.5 h-2.5 rounded-sm" style="background:#ff7d45"></span>Very Low (pLDDT &lt; 50)</span>
                </div>
            </div>
            <div id="ff3dPaneRight" class="relative rounded-xl overflow-hidden border border-slate-700 bg-white" style="height:430px;">
                <div class="absolute top-3 left-3 right-3 z-20 pointer-events-none flex items-start justify-between gap-2 bg-slate-900/90 border border-slate-700/80 rounded-lg px-3 py-2 backdrop-blur-sm shadow-sm">
                    <div>
                        <p id="ff3dRightTitle" class="text-xs font-extrabold text-emerald-400">${is3dp ? 'RIGHT · 3DP Structure-Aware Design (Post-CDS)' : 'RIGHT · 1DP Baseline Design (Post-CDS)'}</p>
                        <p id="ff3dRightSub" class="text-[10px] text-slate-300">${is3dp ? `CAI ${(payload.cai3dp ?? payload.cai1dp).toFixed(3)} · GC ${(payload.gc3dp ?? payload.gc1dp).toFixed(1)}% · 3DP translational pause targets mapped onto identical 3D backbone` : `CAI ${payload.cai1dp.toFixed(3)} · GC ${payload.gc1dp.toFixed(1)}% · 1DP preferred codon substitutions mapped onto identical 3D backbone`}</p>
                    </div>
                    <span id="ff3dRightBadge" class="px-2 py-0.5 rounded-full text-[10px] font-extrabold bg-emerald-950 text-emerald-300 border border-emerald-500/40 shrink-0">${is3dp ? '3DP (Structure-Aware)' : '1DP (Baseline)'}</span>
                </div>
                <div id="ff3dViewerRight" class="w-full h-full relative"><canvas id="ff3dCanvasRight" class="w-full h-full block cursor-grab"></canvas></div>
                <div id="ff3dRightLegend" class="absolute bottom-2.5 left-2.5 right-2.5 z-20 pointer-events-none flex flex-wrap items-center gap-3 bg-slate-900/90 border border-slate-700 rounded-lg px-2.5 py-1.5 text-[10px] text-slate-200 shadow-sm">
                    ${is3dp ? `
                        <span class="inline-flex items-center gap-1"><span class="w-2.5 h-2.5 rounded-full" style="background:#ef4444"></span>3DP Pause Codon Slowdown (Helix-loop boundary)</span>
                        <span class="inline-flex items-center gap-1"><span class="w-2.5 h-2.5 rounded-sm" style="background:#f59e0b"></span>Conserved Pause Boundary</span>
                        <span class="inline-flex items-center gap-1"><span class="w-2.5 h-2.5 rounded-sm" style="background:#10b981"></span>1DP≠3DP Synonymous Shift</span>
                        <span class="inline-flex items-center gap-1"><span class="w-2.5 h-2.5 rounded-sm" style="background:#65cbf3"></span>High-Speed Translation Core</span>
                    ` : `
                        <span class="inline-flex items-center gap-1"><span class="w-2.5 h-2.5 rounded-sm" style="background:#10b981"></span>Preferred Codon Substitution (&gt;42%)</span>
                        <span class="inline-flex items-center gap-1"><span class="w-2.5 h-2.5 rounded-sm" style="background:#0053d6"></span>Synonymous Substitution</span>
                        <span class="inline-flex items-center gap-1"><span class="w-2.5 h-2.5 rounded-sm" style="background:#65cbf3"></span>Unchanged Ribbon Backbone</span>
                    `}
                </div>
            </div>
        </div>
        <div class="mt-3 rounded-xl bg-slate-900/95 border border-slate-800 p-3">
            <div class="flex flex-wrap items-center justify-between gap-2 pb-2 mb-2 border-b border-slate-800 text-[11px]">
                <span class="font-extrabold text-slate-300">Sequence Track (1 ~ ${payload.residues.length} aa) · Click residue to focus 3D structures (Red border = 3DP translation pause target)</span>
                <div id="ff3dInspectorText" class="font-mono text-[11px] text-sky-300"></div>
            </div>
            <div id="ff3dSeqStrip" class="grid gap-1 max-h-36 overflow-y-auto pr-1" style="grid-template-columns: repeat(auto-fill, minmax(42px, 1fr));"></div>
        </div>
    `;

    const inspectorEl = document.getElementById('ff3dInspectorText');
    const seqStripEl = document.getElementById('ff3dSeqStrip');

    function updateInspector(r) {
        if (!inspectorEl || !r) return;
        const origText = r.hasOrigDna ? `${r.codonOrig} (${(r.freqOrig * 100).toFixed(0)}%)` : `${r.codonOrig}`;
        if (is3dp) {
            const statusTag = r.isPause && r.isDiff
                ? '★ 3DP Pause Codon Slowdown'
                : (r.isPause ? 'Conserved Pause Boundary' : (r.isDiff ? '1DP≠3DP Synonymous Shift' : 'Helix Core'));
            inspectorEl.textContent = `Pos ${r.pos} (${r.aa}) | 1DP: ${r.codon1dp} (${(r.freq1dp * 100).toFixed(0)}%) → 3DP: ${r.codon3dp} (${(r.freq3dp * 100).toFixed(0)}%) | pLDDT ${r.plddt.toFixed(1)} | ${statusTag}`;
        } else {
            inspectorEl.textContent = `Pos ${r.pos} (${r.aa}) | Orig: ${origText} → 1DP: ${r.codon1dp} (${(r.freq1dp * 100).toFixed(0)}%) | pLDDT ${r.plddt.toFixed(1)}`;
        }
    }

    function renderSeqStrip() {
        if (!seqStripEl) return;
        seqStripEl.innerHTML = '';
        payload.residues.forEach(r => {
            const cell = document.createElement('button');
            cell.type = 'button';
            const isSel = r.pos === selectedPos;
            const activeCodon = rightMode === '3dp' ? (r.codon3dp || r.codon1dp) : r.codon1dp;
            let bgStyle = 'background:#161f30;color:#94a3b8;';
            if (rightMode === '3dp') {
                if (r.isPause && r.isDiff) bgStyle = 'background:rgba(239,68,68,0.28);color:#fecaca;';
                else if (r.isPause) bgStyle = 'background:rgba(234,179,8,0.22);color:#fef08a;';
                else if (r.isDiff) bgStyle = 'background:rgba(16,185,129,0.22);color:#a7f3d0;';
            } else if (r.isDiffOrig1dp) {
                bgStyle = 'background:rgba(16,185,129,0.22);color:#a7f3d0;';
            }
            const borderStyle = isSel
                ? 'outline:2px solid #ffffff;z-index:2;'
                : (r.isPause ? 'border:1px solid rgba(239,68,68,0.75);' : 'border:1px solid #1e293b;');
            cell.style.cssText = `${bgStyle}${borderStyle}border-radius:5px;padding:2px 1px;font-family:ui-monospace,monospace;font-size:9px;line-height:1.15;text-align:center;cursor:pointer;`;
            cell.innerHTML = `<div style="font-weight:800;">${r.pos}${r.aa}</div><div style="opacity:0.85;">${activeCodon}</div>`;
            cell.addEventListener('click', () => selectResidue(r.pos));
            seqStripEl.appendChild(cell);
        });
    }

    let drawLeft = () => {};
    let drawRight = () => {};
    let viewerLeft = null;
    let viewerRight = null;

    const use3Dmol = Boolean(window.$3Dmol && payload.cartoonPdb);
    if (use3Dmol) {
        const leftHost = document.getElementById('ff3dViewerLeft');
        const rightHost = document.getElementById('ff3dViewerRight');
        leftHost.innerHTML = '';
        rightHost.innerHTML = '';

        viewerLeft = window.$3Dmol.createViewer(leftHost, { backgroundColor: viewerBg, antialias: true });
        viewerRight = window.$3Dmol.createViewer(rightHost, { backgroundColor: viewerBg, antialias: true });

        viewerLeft.addModel(payload.cartoonPdb, 'pdb');
        viewerRight.addModel(payload.cartoonPdb, 'pdb');

        function getAtomCenter(v, pos, fallbackR) {
            const atoms = v.selectedAtoms({ resi: pos, atom: 'CA' });
            if (atoms && atoms.length) {
                return { x: atoms[0].x, y: atoms[0].y, z: atoms[0].z };
            }
            return { x: fallbackR.x, y: fallbackR.y, z: fallbackR.z };
        }

        function styleViewer(v, modeKey) {
            if (!v) return;
            v.removeAllShapes();
            v.removeAllLabels();
            v.setStyle({}, {
                cartoon: {
                    arrows: true,
                    thickness: 0.32,
                    colorfunc: atom => {
                        const r = residueByPos.get(atom.resi);
                        if (!r) return '#0053d6';
                        return modeKey === 'orig' ? ff3dColorLeftInput(r) : ff3dColorRightPost(r, modeKey);
                    }
                }
            });

            // On 3DP panel, add subtle 3D spheres directly on the C-alpha atoms of pause-region codon modifications
            if (modeKey === '3dp') {
                for (const r of payload.residues) {
                    if (r.isPause && r.isDiff) {
                        const centerPt = getAtomCenter(v, r.pos, r);
                        v.addSphere({
                            center: centerPt,
                            radius: 1.25,
                            color: '#ef4444',
                            opacity: 0.90
                        });
                    }
                }
            }

            // Highlight currently selected residue with a callout label + sphere
            const selR = residueByPos.get(selectedPos);
            if (selR) {
                const selCenter = getAtomCenter(v, selR.pos, selR);
                v.addSphere({
                    center: selCenter,
                    radius: 1.95,
                    color: modeKey === 'orig' ? '#0053d6' : (selR.isPause ? '#ef4444' : '#10b981'),
                    opacity: 0.95
                });
                const labelText = modeKey === 'orig'
                    ? `${selR.aa}${selR.pos} · pLDDT ${selR.plddt.toFixed(1)}`
                    : (modeKey === '3dp'
                        ? `${selR.aa}${selR.pos} · 1DP:${selR.codon1dp}→3DP:${selR.codon3dp}`
                        : `${selR.aa}${selR.pos} · 1DP:${selR.codon1dp} (${(selR.freq1dp * 100).toFixed(0)}%)`);
                v.addLabel(labelText, {
                    position: selCenter,
                    backgroundColor: '#0f172a',
                    backgroundOpacity: 0.88,
                    fontColor: '#ffffff',
                    fontSize: 11,
                    borderThickness: 1,
                    borderColor: modeKey === 'orig' ? '#38bdf8' : '#10b981'
                });
            }

            v.setClickable({}, true, atom => {
                if (atom && atom.resi) selectResidue(atom.resi);
            });
            v.render();
        }

        function resetBothToCanonicalView() {
            if (!viewerLeft || !viewerRight) return;
            styleViewer(viewerLeft, leftMode);
            styleViewer(viewerRight, rightMode);
            viewerLeft.zoomTo();
            viewerLeft.rotate(-68, 'y');
            viewerLeft.rotate(14, 'x');
            viewerLeft.zoom(1.06);
            viewerLeft.render();
            viewerRight.setView(viewerLeft.getView());
            viewerRight.render();
        }

        drawLeft = () => styleViewer(viewerLeft, leftMode);
        drawRight = () => styleViewer(viewerRight, rightMode);
        resetBothToCanonicalView();

        // Optional camera sync when syncCamera is explicitly enabled by user
        let syncingView = false;
        const syncLeftToRight = () => {
            if (!syncCamera || syncingView || !viewerLeft || !viewerRight) return;
            syncingView = true;
            viewerRight.setView(viewerLeft.getView());
            viewerRight.render();
            syncingView = false;
        };
        const syncRightToLeft = () => {
            if (!syncCamera || syncingView || !viewerLeft || !viewerRight) return;
            syncingView = true;
            viewerLeft.setView(viewerRight.getView());
            viewerLeft.render();
            syncingView = false;
        };
        leftHost.addEventListener('mousemove', e => { if (e.buttons) syncLeftToRight(); });
        leftHost.addEventListener('wheel', () => setTimeout(syncLeftToRight, 10), { passive: true });
        rightHost.addEventListener('mousemove', e => { if (e.buttons) syncRightToLeft(); });
        rightHost.addEventListener('wheel', () => setTimeout(syncRightToLeft, 10), { passive: true });
    } else {
        let camLeft = { rx: 0.22, ry: -1.15, zoom: 1.0, panX: 0, panY: 0 };
        let camRight = { rx: 0.22, ry: -1.15, zoom: 1.0, panX: 0, panY: 0 };
        drawLeft = attachInteractive3dCanvas(
            document.getElementById('ff3dCanvasLeft'),
            payload.residues,
            () => leftMode,
            () => camLeft,
            next => {
                camLeft = next;
                if (syncCamera) camRight = { ...next };
                drawLeft();
                if (syncCamera) drawRight();
            },
            () => selectedPos,
            pos => selectResidue(pos),
            r => { if (r) updateInspector(r); }
        );
        drawRight = attachInteractive3dCanvas(
            document.getElementById('ff3dCanvasRight'),
            payload.residues,
            () => rightMode,
            () => camRight,
            next => {
                camRight = next;
                if (syncCamera) camLeft = { ...next };
                drawRight();
                if (syncCamera) drawLeft();
            },
            () => selectedPos,
            pos => selectResidue(pos),
            r => { if (r) updateInspector(r); }
        );
    }

    function selectResidue(pos) {
        selectedPos = pos;
        const found = residueByPos.get(pos);
        if (found) updateInspector(found);
        renderSeqStrip();
        drawLeft();
        drawRight();
    }

    function updatePaneLabels() {
        const lTitle = document.getElementById('ff3dLeftTitle');
        const lSub = document.getElementById('ff3dLeftSub');
        const lBadge = document.getElementById('ff3dLeftBadge');
        const rTitle = document.getElementById('ff3dRightTitle');
        const rSub = document.getElementById('ff3dRightSub');
        const rBadge = document.getElementById('ff3dRightBadge');
        const rLegend = document.getElementById('ff3dRightLegend');

        if (leftMode === 'orig') {
            lTitle.textContent = 'LEFT · Original Structure (Input Structure)';
            lSub.textContent = 'Input coordinates · 3D cartoon ribbon (pLDDT confidence colors)';
            lBadge.textContent = 'Input coordinates';
        } else {
            lTitle.textContent = 'LEFT · 1DP Baseline Design (Baseline 1DP)';
            lSub.textContent = `CAI ${payload.cai1dp.toFixed(3)} · GC ${payload.gc1dp.toFixed(1)}% · Uniform preferred codons across all domains`;
            lBadge.textContent = '1DP Baseline';
        }

        if (rightMode === '3dp') {
            rTitle.textContent = 'RIGHT · 3DP Structure-Aware Design (Post-CDS)';
            rSub.textContent = `CAI ${(payload.cai3dp ?? payload.cai1dp).toFixed(3)} · GC ${(payload.gc3dp ?? payload.gc1dp).toFixed(1)}% · 3DP translational pause targets mapped onto identical 3D backbone`;
            rBadge.textContent = '3DP (Structure-Aware)';
            rLegend.innerHTML = `
                <span class="inline-flex items-center gap-1"><span class="w-2.5 h-2.5 rounded-full" style="background:#ef4444"></span>3DP Pause Codon Slowdown (Helix-loop boundary)</span>
                <span class="inline-flex items-center gap-1"><span class="w-2.5 h-2.5 rounded-sm" style="background:#f59e0b"></span>Conserved Pause Boundary</span>
                <span class="inline-flex items-center gap-1"><span class="w-2.5 h-2.5 rounded-sm" style="background:#10b981"></span>1DP≠3DP Synonymous Shift</span>
                <span class="inline-flex items-center gap-1"><span class="w-2.5 h-2.5 rounded-sm" style="background:#65cbf3"></span>High-Speed Translation Core</span>
            `;
        } else {
            rTitle.textContent = 'RIGHT · 1DP Baseline Design (Post-CDS)';
            rSub.textContent = `CAI ${payload.cai1dp.toFixed(3)} · GC ${payload.gc1dp.toFixed(1)}% · 1DP preferred codon substitutions mapped onto identical 3D backbone`;
            rBadge.textContent = '1DP (Baseline)';
            rLegend.innerHTML = `
                <span class="inline-flex items-center gap-1"><span class="w-2.5 h-2.5 rounded-sm" style="background:#10b981"></span>Preferred Codon Substitution (&gt;42%)</span>
                <span class="inline-flex items-center gap-1"><span class="w-2.5 h-2.5 rounded-sm" style="background:#0053d6"></span>Synonymous Substitution</span>
                <span class="inline-flex items-center gap-1"><span class="w-2.5 h-2.5 rounded-sm" style="background:#65cbf3"></span>Unchanged Ribbon Backbone</span>
            `;
        }
    }

    container.querySelectorAll('[data-3d-pair]').forEach(btn => {
        btn.addEventListener('click', () => {
            const pair = btn.getAttribute('data-3d-pair');
            container.querySelectorAll('[data-3d-pair]').forEach(b => {
                b.className = b === btn
                    ? 'px-3 py-1.5 rounded-lg text-[11px] font-extrabold bg-emerald-600 text-white shadow-sm'
                    : 'px-3 py-1.5 rounded-lg text-[11px] font-bold bg-slate-800 text-slate-300 hover:bg-slate-700';
            });
            if (pair === 'orig_3dp') {
                leftMode = 'orig'; rightMode = '3dp';
            } else if (pair === 'orig_1dp') {
                leftMode = 'orig'; rightMode = '1dp';
            } else {
                leftMode = '1dp'; rightMode = '3dp';
            }
            updatePaneLabels();
            renderSeqStrip();
            drawLeft();
            drawRight();
        });
    });

    const bgBtn = document.getElementById('ff3dBgBtn');
    bgBtn?.addEventListener('click', () => {
        viewerBg = viewerBg === '#ffffff' ? '#090d16' : '#ffffff';
        bgBtn.textContent = viewerBg === '#ffffff' ? 'Background: White (AlphaFold)' : 'Background: Dark mode';
        const pL = document.getElementById('ff3dPaneLeft');
        const pR = document.getElementById('ff3dPaneRight');
        if (pL) pL.style.backgroundColor = viewerBg;
        if (pR) pR.style.backgroundColor = viewerBg;
        if (viewerLeft && viewerRight) {
            viewerLeft.setBackgroundColor(viewerBg);
            viewerRight.setBackgroundColor(viewerBg);
            viewerLeft.render();
            viewerRight.render();
        }
    });

    const syncBtn = document.getElementById('ff3dSyncBtn');
    syncBtn?.addEventListener('click', () => {
        syncCamera = !syncCamera;
        syncBtn.textContent = `Camera Sync: ${syncCamera ? 'On (Synchronized)' : 'Off (Independent)'}`;
        syncBtn.className = syncCamera
            ? 'px-2.5 py-1.5 rounded-lg bg-sky-900/70 hover:bg-sky-800 text-sky-200 text-[10px] font-bold border border-sky-500/40'
            : 'px-2.5 py-1.5 rounded-lg bg-slate-800 hover:bg-slate-700 text-slate-200 text-[10px] font-bold border border-slate-700';
    });

    document.getElementById('ff3dResetBtn')?.addEventListener('click', () => {
        if (viewerLeft && viewerRight) {
            viewerLeft.zoomTo();
            viewerLeft.rotate(-68, 'y');
            viewerLeft.rotate(14, 'x');
            viewerLeft.zoom(1.06);
            viewerLeft.render();
            viewerRight.setView(viewerLeft.getView());
            viewerRight.render();
        } else {
            drawLeft();
            drawRight();
        }
    });

    const spinBtn = document.getElementById('ff3dSpinBtn');
    spinBtn?.addEventListener('click', () => {
        spin = !spin;
        spinBtn.textContent = `Auto Spin: ${spin ? 'On' : 'Off'}`;
        if (viewerLeft && viewerRight) {
            viewerLeft.spin(spin ? 'y' : false, 0.5);
            viewerRight.spin(spin ? 'y' : false, 0.5);
        }
    });

    const initFound = residueByPos.get(selectedPos);
    if (initFound) updateInspector(initFound);
    renderSeqStrip();
    if (!use3Dmol) {
        drawLeft();
        drawRight();
    }
}

function renderResults() {
    renderScientificReview(state.results || {});
    const res = state.results;
    if (!res) return;
    const primary = getPrimaryResult(res);
    renderExecutionEvidence(res);

    elements.emptyState.classList.add('hidden');
    elements.resultsContainer.classList.remove('hidden');
    setWorkspaceState('review');
    if (elements.resultContextSummary) {
        const host = formatHostProfile(getResultHostProfile(res));
        elements.resultContextSummary.textContent = `${host} · Review the computational checks before synthesis or experimental use.`;
    }
    if (elements.designContractSummary) {
        const contract = res.design_contract;
        if (contract?.engine_id === 'dp_v2_1_1') {
            const axes = (contract.scientific_axes || [])
                .map(axis => `<span class="rounded-lg border border-teal-200 bg-white px-2 py-1 font-bold text-teal-800 dark:border-teal-800 dark:bg-slate-900 dark:text-teal-200">${escapeHtml(axis.id)} · ${escapeHtml(axis.evidence_class)}</span>`)
                .join('');
            const localGuard = contract.local_composition_guard || {};
            const activeBand = Array.isArray(localGuard.initiation_gc_active_count_band)
                ? localGuard.initiation_gc_active_count_band.join('–')
                : 'not reported';
            elements.designContractSummary.innerHTML = `<div class="font-extrabold text-teal-900 dark:text-teal-100">DP v2.1.1 ${escapeHtml(contract.engine_version)} · development candidate</div><div class="mt-3 flex flex-wrap gap-2">${axes}</div><p class="mt-3 text-slate-600 dark:text-slate-300">Local guard: active 5′ GC count layer ${escapeHtml(activeBand)} · homopolymer ceiling ${escapeHtml(localGuard.homopolymer_max_run ?? 'not reported')} nt. RNA folding remains an independently evaluated metric.</p>`;
            elements.designContractSummary.classList.remove('hidden');
        } else {
            elements.designContractSummary.innerHTML = '';
            elements.designContractSummary.classList.add('hidden');
        }
    }

    if (res.construct_id) {
        elements.constructIdDisplay.textContent = res.construct_id;
        elements.constructIdRow.classList.remove('hidden');
    } else {
        elements.constructIdDisplay.textContent = '';
        elements.constructIdRow.classList.add('hidden');
    }

    // Metrics
    elements.caiValue.textContent = primary.metrics.cai.toFixed(3);
    // GC Value updated below via manual calculation
    elements.polyaValue.textContent = primary.metrics.polya_signals === 0 ? '0 (Clean)' : primary.metrics.polya_signals;
    if (elements.hostProfileValue) {
        elements.hostProfileValue.textContent = formatHostProfile(getResultHostProfile(res));
    }

    // Metrics Comparison Table
    const isProteinInput = isProteinInputResult(res);
    elements.origLen.textContent = `${res.original_length || state.sequence.length} ${isProteinInput ? 'aa' : 'bp'}`;
    elements.optLen.textContent = `${primary.metrics.length || primary.optimized_sequence.length} bp`;

    // Variables for calculations
    const origSeq = state.sequence;
    const hasOriginalSequence = Boolean(origSeq);
    const optSeq = primary.optimized_sequence;
    const calculatedGC = calculateGC(optSeq);
    const oGC = isProteinInput || !hasOriginalSequence ? null : calculateGC(origSeq);

    elements.origGC.textContent = isProteinInput ? 'N/A' : (hasOriginalSequence ? `${oGC}%` : 'Not recorded');
    elements.optGCComp.textContent = `${calculatedGC.toFixed(1)}%`;
    elements.gcValue.textContent = `${calculatedGC.toFixed(1)}%`;
    const gcTarget = getResultGcTarget(res, primary);
    const baseReportData = buildResultsReportModel(res, primary, gcTarget);
    const cand1dpSeq = ((Array.isArray(res.candidates) ? res.candidates : (res.report_candidates || [])).find(c => c.id === 'baseline_1dp')?.dna_sequence || '');
    const refSeqForAccounting = (!isProteinInput && origSeq)
        ? origSeq
        : (cand1dpSeq || synthesizeBaselineReferenceDna(optSeq, res));
    const accounting = buildCodonAccounting(refSeqForAccounting, optSeq, res, gcTarget);
    const structureViewerPayload = buildStructureViewerPayload(res, primary);
    state.reportData = { ...baseReportData, accounting, codons: accounting.codons, structure_viewer: structureViewerPayload };
    renderCodonAccounting(accounting);
    renderCodonTracks(accounting);
    renderStructure3dViewer(res, primary);
    elements.gcTargetRange.textContent = `Target: ${gcTarget.min.toFixed(1)}–${gcTarget.max.toFixed(1)}%`;
    const originalEvaluation = res.acceptance_evaluation?.original;
    const originalCai = originalEvaluation?.criteria?.find(row => row.criterion === 'cai')?.observed;
    elements.origCAI.textContent = originalCai == null ? 'Unavailable' : Number(originalCai).toFixed(3);
    elements.optCAIComp.textContent = primary.metrics.cai.toFixed(3);

    const mutationRow = elements.mutationRate?.closest('tr, .metric-row');
    if (isProteinInput) {
        elements.mutationRate.textContent = 'N/A';
        if (mutationRow) mutationRow.classList.add('hidden');
    } else if (!hasOriginalSequence) {
        if (mutationRow) mutationRow.classList.remove('hidden');
        elements.mutationRate.textContent = 'Not recorded';
    } else {
        if (mutationRow) mutationRow.classList.remove('hidden');
        // Calculate Mutation Rate
        let diffCount = 0;
        const compareLen = Math.min(origSeq.length, optSeq.length);
        for (let i = 0; i < compareLen; i++) {
            if (origSeq[i] !== optSeq[i]) diffCount++;
        }
        diffCount += Math.abs(origSeq.length - optSeq.length);
        const mRate = origSeq.length > 0 ? ((diffCount / Math.max(origSeq.length, 1)) * 100).toFixed(1) : 0;
        elements.mutationRate.textContent = `${mRate}% (${diffCount} bp)`;
    }
    if (elements.aaIdentity) elements.aaIdentity.textContent = primary.aaPreserved;
    renderDesignReview(res);

    // Sequence
    elements.optimizedSequence.textContent = formatSequence(primary.optimized_sequence);

    // Render GC Graph
    renderGCGraph(primary.optimized_sequence, getResultHostProfile(res));
    renderCandidateComparison(res);
    renderCustomRestrictionResults(res);
    renderMfeWarning(res);
    renderResultsReport(res, primary, gcTarget, state.reportData);
    renderComparisonDashboard(res);

    // PolyA color coding
    const polyaCount = primary.metrics.polya_signals;
    if (polyaCount === 0) {
        elements.polyaValue.textContent = '0 (Clean)';
        elements.polyaValue.className = 'text-xl font-black text-emerald-600 dark:text-emerald-400 leading-none';
    } else {
        elements.polyaValue.textContent = `${polyaCount} ⚠️`;
        elements.polyaValue.className = 'text-xl font-black text-amber-600 dark:text-amber-400 leading-none';
    }

    // Validation Badges (registry-driven)
    elements.validationStatus.classList.remove('hidden');
    renderValidationChecks(validationRegistry, resultChecksFromResponse(res));

    // Structure prediction links
    updateStructureLinks();

    // JSON Details
    elements.jsonDetails.textContent = JSON.stringify(res, null, 2);
}

// Report only response-backed evidence, never a fixed successful pipeline.
function renderExecutionEvidence(res) {
    const grid = document.getElementById('traceStagesGrid');
    if (!grid) return;
    const contract = res.design_contract;
    const profile = res.sop_profile;
    const stages = [
        ['Design response', 'Received', 'A response was returned; this is not biological validation.'],
        ['Sequence policy', profile?.profile_hash ? 'Reported' : 'Not reported', profile?.profile_hash ? `${profile.profile_name || profile.profile_id || 'Profile'} · ${profile.profile_version || 'version not reported'}. See response JSON for the profile hash.` : 'Applied profile provenance is unavailable in this response.'],
        ['Design engine', contract?.engine_id ? 'Reported' : 'Not reported', contract?.engine_id ? `${contract.engine_id} · ${contract.engine_version || 'version not reported'}` : 'Engine provenance is unavailable in this response.'],
        ['Sequence checks', 'Inspect results', 'See individual sequence checks; execution does not imply that all checks passed.'],
        ['Policy disposition', 'Inspect results', 'See the QC decision matrix for the actual candidate disposition.'],
        ['Platform persistence / RNA analysis', 'Not verified here', 'This panel does not verify storage, agent execution, or RNA analysis.']
    ];
    grid.innerHTML = stages.map(([title, status, detail]) => `<div class="p-3 rounded-xl border border-slate-200 dark:border-slate-800"><div class="flex justify-between text-xs font-bold"><span>${escapeHtml(title)}</span><span>${escapeHtml(status)}</span></div><p class="text-[10px] text-slate-500 mt-1">${escapeHtml(detail)}</p></div>`).join('');
}

const ALIGNMENT_PAGE_SIZE = 60;

function comparisonNumber(value, minimum, maximum) {
    return typeof value === 'number' && Number.isFinite(value)
        && value >= minimum && value <= maximum ? value : null;
}

function comparisonIdentity(metrics) {
    const identity = comparisonNumber(metrics.aa_identity, 0, 1);
    if (identity === null) return { text: 'Not evaluated', status: 'unknown' };
    const passed = identity === 1;
    return { text: `${(identity * 100).toFixed(2)}% ${passed ? 'Passed' : 'Failed'}`, status: passed ? 'pass' : 'fail' };
}

function comparisonTypeIIS(metrics) {
    const clean = metrics.type_iis_clean;
    const count = metrics.type_iis_site_count;
    const validCount = Number.isInteger(count) && count >= 0;
    if (clean === false) {
        return { text: validCount && count > 0 ? `${count} site(s) — Failed` : 'Failed', status: 'fail' };
    }
    if (clean === true && (count === undefined || (validCount && count === 0))) {
        return { text: 'Clean', status: 'pass' };
    }
    return { text: 'Not evaluated', status: 'unknown' };
}

function comparisonStatus(rule, ml) {
    if (rule.status === 'fail' || ml.status === 'fail') return 'Failed';
    return rule.status === 'pass' && ml.status === 'pass' ? 'Passed' : 'Not evaluated';
}

function renderComparisonDashboard(res) {
    const comparison = res?.comparison;
    if (!comparison || res.mode !== 'dual_compare') {
        elements.comparisonDashboard.classList.add('hidden');
        return;
    }
    elements.comparisonDashboard.classList.remove('hidden');
    elements.comparisonNotice.textContent = res.preview_notice || 'Experimental comparison preview.';
    const rule = comparison.rule?.metrics || {};
    const ml = comparison.ml?.metrics || {};
    const ruleIdentity = comparisonIdentity(rule);
    const mlIdentity = comparisonIdentity(ml);
    const ruleTypeIIS = comparisonTypeIIS(rule);
    const mlTypeIIS = comparisonTypeIIS(ml);
    const numericRow = (label, key, max, digits, suffix, deltaSuffix) => {
        const left = comparisonNumber(rule[key], 0, max);
        const right = comparisonNumber(ml[key], 0, max);
        return [label,
            left === null ? 'Not evaluated' : `${left.toFixed(digits)}${suffix}`,
            right === null ? 'Not evaluated' : `${right.toFixed(digits)}${suffix}`,
            left === null || right === null ? 'Not evaluated' : `${(right - left).toFixed(digits)}${deltaSuffix}`];
    };
    const percentText = (value, suffix) => {
        const number = comparisonNumber(value, 0, 100);
        return number === null ? 'Not evaluated' : `${number.toFixed(1)}%${suffix}`;
    };
    const metricRows = [
        ['AA Translation Identity', ruleIdentity.text, mlIdentity.text, comparisonStatus(ruleIdentity, mlIdentity)],
        numericRow('GC Content', 'gc_percent', 100, 1, '%', ' pp'),
        numericRow('CAI Index', 'cai', 1, 3, '', ''),
        ['Type IIS Clearance', ruleTypeIIS.text, mlTypeIIS.text, comparisonStatus(ruleTypeIIS, mlTypeIIS)],
        ['Codon Concordance', '—', '—', percentText(comparison.codon_concordance_percent, ' match')],
        ['NT Identity', '—', '—', percentText(comparison.nt_identity_percent, '')]
    ];
    elements.comparisonMatrixBody.innerHTML = metricRows.map(row => `<tr>${row.map((cell, index) => `<td class="p-3 ${index === 0 ? 'font-bold text-slate-700 dark:text-slate-200' : 'text-slate-600 dark:text-slate-300'}">${escapeHtml(String(cell))}</td>`).join('')}</tr>`).join('');

    const provenance = comparison.provenance || {};
    const badge = (label, status) => {
        const verified = status === 'verified';
        const classes = verified ? 'bg-emerald-50 text-emerald-700 border-emerald-200' : 'bg-slate-100 text-slate-600 border-slate-200';
        return `<span class="px-2 py-1 rounded-lg border text-[10px] font-bold ${classes}">${escapeHtml(label)}: ${escapeHtml(status || 'unavailable')}</span>`;
    };
    elements.provenanceBadges.innerHTML = [
        badge('Canonical DB', provenance.db_save_status),
        badge('Audit immutability', provenance.audit_status),
        badge('Evaluation leakage', provenance.leakage_check_status),
        provenance.canonical_sha256 ? `<button class="px-2 py-1 rounded-lg border border-indigo-200 bg-indigo-50 text-indigo-700 text-[10px] font-mono" title="${escapeHtml(provenance.canonical_sha256)}" onclick="navigator.clipboard.writeText('${escapeHtml(provenance.canonical_sha256)}')">SHA-256 · Copy</button>` : ''
    ].join('');
    renderCodonAlignment(comparison.alignment || []);
}

function renderCodonAlignment(alignment) {
    const maxPage = Math.max(0, Math.ceil(alignment.length / ALIGNMENT_PAGE_SIZE) - 1);
    state.alignmentPage = Math.min(Math.max(state.alignmentPage, 0), maxPage);
    const start = state.alignmentPage * ALIGNMENT_PAGE_SIZE;
    const points = alignment.slice(start, start + ALIGNMENT_PAGE_SIZE);
    elements.codonAlignmentViewer.style.setProperty('--codon-count', points.length);
    const row = (label, formatter, className = '') => [
        `<div class="codon-label">${label}</div>`,
        ...points.map(point => `<div class="codon-chip ${className && formatter(point, true) ? className : ''}" title="${escapeHtml(codonTooltip(point))}">${escapeHtml(String(formatter(point, false)))}</div>`)
    ].join('');
    elements.codonAlignmentViewer.innerHTML = [
        row('AA', point => point.amino_acid),
        row('Rule', point => point.rule_codon),
        row('ML', (point, classProbe) => classProbe ? point.is_different : point.ml_codon, 'ml-different'),
        row('Diff', point => point.is_different ? '▲' : '·', 'diff-marker')
    ].join('');
    elements.alignmentRange.textContent = alignment.length
        ? `Codons ${start + 1}–${Math.min(start + points.length, alignment.length)} of ${alignment.length}`
        : 'No alignment data';
    elements.alignmentPrev.disabled = state.alignmentPage === 0;
    elements.alignmentNext.disabled = state.alignmentPage >= maxPage;
}

function codonTooltip(point) {
    return `AA ${point.amino_acid} · position ${point.position}\nRule ${point.rule_codon}: GC ${point.rule_gc_bases}/3, host frequency ${(Number(point.rule_frequency) * 100).toFixed(2)}%\nML ${point.ml_codon}: GC ${point.ml_gc_bases}/3, host frequency ${(Number(point.ml_frequency) * 100).toFixed(2)}%`;
}

function changeAlignmentPage(delta) {
    const alignment = state.results?.comparison?.alignment || [];
    if (!alignment.length) return;
    state.alignmentPage += delta;
    renderCodonAlignment(alignment);
}

function getAcceptanceCriteriaPayload() {
    return {
        cai: { mode: elements.criterionCaiMode.value, minimum: 0.8 },
        overall_gc: { mode: elements.criterionGcMode.value },
        local_gc: { mode: elements.criterionLocalGcMode.value, minimum: 30, maximum: 70, window_size: 60 },
        type_iis: { mode: elements.criterionTypeIisMode.value, enzymes: ['BsaI', 'BsmBI/Esp3I', 'SapI'], custom_sites: state.customRestrictionSites },
        repeats: { mode: elements.criterionRepeatsMode.value, minimum_length: 18, maximum_count: 0 },
        homopolymers: { mode: elements.criterionHomopolymerMode.value, maximum_length: 8 },
        forbidden_motifs: { mode: elements.criterionMotifsMode.value, motifs: ['AATAAA', 'GTAAGT', 'ATTTA'] }
    };
}

function renderDesignReview(res) {
    const decision = res.automated_decision || 'Unavailable';
    const summary = res.decision_summary || {};
    elements.automatedDecisionValue.textContent = decision.replace('_', ' ');
    elements.automatedDecisionValue.setAttribute('aria-label', `Automated decision: ${decision}`);
    elements.automatedDecisionSummary.textContent = res.decision_summary
        ? `${summary.required_failure_count || 0} required failures · ${summary.preferred_warning_count || 0} preferred warnings. ${summary.explanation || ''}`
        : 'Acceptance criteria were not returned by the API.';
    const rows = res.qc_decision_matrix || [];
    if (!rows.length) {
        elements.qcDecisionMatrix.classList.add('hidden');
    } else {
        elements.qcDecisionMatrix.classList.remove('hidden');
        elements.qcDecisionMatrixBody.innerHTML = `<table class="w-full text-xs"><thead class="bg-slate-50"><tr><th class="px-3 py-2 text-left">Criterion</th><th class="px-3 py-2 text-left">Mode</th><th class="px-3 py-2 text-left">Observed</th><th class="px-3 py-2 text-left">Result</th></tr></thead><tbody class="divide-y divide-slate-100">${rows.map(row => `<tr><td class="px-3 py-2">${escapeHtml(row.criterion)}</td><td class="px-3 py-2">${escapeHtml(row.mode)}</td><td class="px-3 py-2">${escapeHtml(String(row.observed))}</td><td class="px-3 py-2 font-bold">${row.result === 'PASS' ? '✓ PASS' : row.result === 'FAIL' ? '✕ FAIL' : row.result === 'WARN' ? '⚠ WARN' : '— IGNORED'}</td></tr>`).join('')}</tbody></table>`;
    }
    const disposition = res.reviewer_disposition;
    elements.reviewerDispositionStatus.textContent = disposition
        ? `${disposition.final_state}: ${disposition.reason || 'No written reason'}`
        : '';
}

function saveReviewerDisposition() {
    if (!state.results) return;
    const disposition = elements.reviewerDisposition.value;
    if (!disposition) {
        showToast('Choose a reviewer disposition first', 'error');
        return;
    }
    const reason = elements.reviewerReason.value.trim();
    if (state.results.automated_decision === 'FAIL' && (disposition === 'accept' || disposition === 'accept_with_exception') && !reason) {
        showToast('A reason is required to accept an automated FAIL', 'error');
        return;
    }
    const finalState = disposition === 'accept' || disposition === 'accept_with_exception'
        ? (state.results.automated_decision === 'FAIL' ? 'MANUALLY_ACCEPTED' : 'ACCEPTED')
        : disposition === 'return_for_redesign' ? 'RETURNED_FOR_REDESIGN' : 'REJECTED';
    state.results.reviewer_disposition = { disposition, reason: reason || null, final_state: finalState, automated_decision: state.results.automated_decision, timestamp: new Date().toISOString() };
    renderDesignReview(state.results);
    addToHistory(state.sequence, state.results);
    showToast('Reviewer disposition saved locally', 'success');
}

function parseCustomRestrictionSites(rawValue) {
    const raw = (rawValue || '').trim();
    if (!raw) return [];

    const entries = raw.split(/[\n,]+/).map(item => item.trim()).filter(Boolean);
    const seenNames = new Set();
    const seenSequences = new Set();

    return entries.map((entry, index) => {
        const parts = entry.split(':');
        const hasName = parts.length > 1;
        const name = hasName ? parts[0].trim() : `Site ${index + 1}`;
        const sequence = (hasName ? parts.slice(1).join(':') : parts[0]).trim().toUpperCase();

        if (!name) {
            throw new Error('Custom restriction site name is required');
        }
        if (!/^[ACGT]+$/.test(sequence)) {
            throw new Error(`Invalid custom restriction site sequence: ${sequence || '(empty)'}`);
        }
        if (sequence.length < 4 || sequence.length > 12) {
            throw new Error(`${name} must be 4-12 bp`);
        }
        if (seenNames.has(name)) {
            throw new Error(`Duplicate custom restriction site name: ${name}`);
        }
        if (seenSequences.has(sequence)) {
            throw new Error(`Duplicate custom restriction site sequence: ${sequence}`);
        }

        seenNames.add(name);
        seenSequences.add(sequence);
        return { name, sequence, scan_rc: true };
    });
}

function mergeRestrictionSitePresets(customSites, selectedEnzymes) {
    const merged = [...customSites];
    const names = new Set(customSites.map(site => site.name.toLowerCase()));
    selectedEnzymes.forEach(name => {
        if (!Object.hasOwn(TYPE_IIS_PRESETS, name) || names.has(name.toLowerCase())) return;
        merged.push({ name, sequence: TYPE_IIS_PRESETS[name], scan_rc: true });
        names.add(name.toLowerCase());
    });
    return merged;
}

function getPrimaryResult(res) {
    const candidate = res.recommended_candidate || (Array.isArray(res.candidates) ? res.candidates[0] : null);
    if (!candidate) {
        return {
            optimized_sequence: res.optimized_sequence,
            metrics: res.metrics,
            validation: res.validation,
            aaPreserved: '✅ 100%'
        };
    }

    const gcWinMin = candidate.gc_window_min != null ? Number(candidate.gc_window_min) : 40.0;
    const gcWinMax = candidate.gc_window_max != null ? Number(candidate.gc_window_max) : 55.0;
    return {
        optimized_sequence: candidate.dna_sequence,
        metrics: {
            cai: Number(candidate.cai || 0),
            gc_percent: Number(candidate.gc_percent || 0),
            polya_signals: Number(candidate.polya_signals || 0),
            length: candidate.dna_sequence ? candidate.dna_sequence.length : 0
        },
        validation: {
            polya: 'PASS',
            moclo: 'UNCHECKED',
            gc: candidate.gc_percent >= gcWinMin && candidate.gc_percent <= gcWinMax ? 'PASS' : 'WARNING'
        },
        gc_window_min: gcWinMin,
        gc_window_max: gcWinMax,
        aaPreserved: res.constraint_report && res.constraint_report.aa_identity === 1.0
            ? '✅ 100%' : '⚠️ Review'
    };
}

// Same wording as the profile-selection cards above (Optimization Settings),
// reused here so the comparison table explains each candidate strategy
// instead of showing a bare name.
const CANDIDATE_DESCRIPTIONS = {
    feasibility_best: 'Dynamic-programming candidate with maximum CAI inside the requested GC range when feasible.',
    gc_target: 'Profile-engine comparison that prioritizes proximity to the active host GC band.',
    high_cai: 'Profile-engine comparison that favors high codon adaptation index.'
};

function renderCandidateComparison(res) {
    if (!elements.candidateComparisonContainer || !elements.candidateComparisonBody) return;
    if (!Array.isArray(res.candidates) || res.candidates.length === 0) {
        elements.candidateComparisonContainer.classList.add('hidden');
        elements.candidateComparisonBody.innerHTML = '';
        return;
    }

    let structureSummary = document.getElementById('structureComparisonSummary');
    if (!structureSummary) {
        structureSummary = document.createElement('p');
        structureSummary.id = 'structureComparisonSummary';
        structureSummary.className = 'mt-3 text-xs text-slate-600 dark:text-slate-300';
        elements.candidateComparisonContainer.appendChild(structureSummary);
    }
    structureSummary.hidden = !res.structure_comparison;
    if (res.structure_comparison) {
        const comparison = res.structure_comparison;
        structureSummary.textContent = `3DP structure-aware design: ΔCAI ${Number(comparison.delta_1d_cai).toFixed(4)} · 3D codon fitness gain +${Number(comparison.delta_3d_gain).toFixed(4)} · Pause-modulated codons: ${comparison.total_edits}`;
    }
    const recommendedId = res.recommended_candidate ? res.recommended_candidate.id : res.candidates[0].id;
    const visibleCandidates = res.design_mode === 'paired_3dp'
        ? res.candidates.filter(c => c.id === 'structure_3dp' || c.id === recommendedId)
        : res.candidates.filter(c => c.id === recommendedId);
    const rowsToRender = visibleCandidates.length ? visibleCandidates : [res.candidates[0]];
    elements.candidateComparisonBody.innerHTML = rowsToRender.map(candidate => {
        const isRecommended = candidate.id === recommendedId;
        const status = candidate.automated_decision || (candidate.validator_status === 'pass' ? 'PASS' : 'REVIEW');
        const aaStatus = candidate.validator_status === 'pass' ? '100%' : 'Review';
        const description = CANDIDATE_DESCRIPTIONS[candidate.id];
        const nameCell = description
            ? `${escapeHtml(candidate.label || candidate.id)}<span class="tooltip-container"><span class="info-icon">i</span><span class="tooltip-text">${escapeHtml(description)}</span></span>`
            : escapeHtml(candidate.label || candidate.id);
        return `
            <tr class="${isRecommended ? 'bg-emerald-50 dark:bg-emerald-900/20' : 'bg-white dark:bg-slate-900'}">
                <td class="px-3 py-2 font-bold text-slate-800 dark:text-slate-100">${nameCell}${isRecommended ? (res.design_mode === 'paired_3dp' ? ' (3DP Selected)' : ' (1DP Selected)') : ''}</td>
                <td class="px-3 py-2">${Number(candidate.cai || 0).toFixed(3)}</td>
                <td class="px-3 py-2">${Number(candidate.gc_percent || 0).toFixed(1)}</td>
                <td class="px-3 py-2">${Number(candidate.gc_window_min || 0).toFixed(1)}-${Number(candidate.gc_window_max || 0).toFixed(1)}</td>
                <td class="px-3 py-2">${aaStatus}</td>
                <td class="px-3 py-2">${candidate.internal_stop_count || 0}</td>
                <td class="px-3 py-2">${candidate.repeat_count || 0}</td>
                <td class="px-3 py-2 font-bold">${status}</td>
            </tr>
        `;
    }).join('');
    elements.candidateComparisonContainer.classList.remove('hidden');
}

function getResultGcTarget(res, primary) {
    const metrics = res.metrics || {};
    const fallback = getGcRange(getResultHostProfile(res));
    return {
        min: Number(metrics.requested_gc_min_percent ?? primary.gc_window_min ?? fallback.gc_min),
        max: Number(metrics.requested_gc_max_percent ?? primary.gc_window_max ?? fallback.gc_max)
    };
}

function getMfeStatus(res) {
    const metrics = res.metrics || {};
    const status = metrics.mfe_status || 'not_computed';
    return {
        status,
        reason: metrics.mfe_status_reason ?? (status === 'computed' ? null : 'status unavailable')
    };
}

function renderMfeWarning(res) {
    if (!elements.mfeWarningBanner) return;
    const mfe = getMfeStatus(res);
    if (mfe.status === 'computed') {
        elements.mfeWarningBanner.textContent = '';
        elements.mfeWarningBanner.classList.add('hidden');
        return;
    }
    elements.mfeWarningBanner.textContent = res.scientific_review
        ? `Full-CDS RNA/MFE analysis was not performed (${mfe.reason}). The separate initiation-window observation is reported below.`
        : `RNA secondary-structure/MFE analysis was not performed (${mfe.reason}).`;
    elements.mfeWarningBanner.classList.remove('hidden');
}

const REPORT_CHECK_LABELS = Object.freeze({
    cai: 'Codon Adaptation Index (CAI)',
    overall_gc: 'Overall GC content',
    local_gc: 'Local GC window',
    type_iis: 'Type IIS sites',
    repeats: 'Direct repeats',
    homopolymers: 'Homopolymers',
    forbidden_motifs: 'Forbidden motifs',
    mfe: 'RNA secondary structure / MFE',
});

const REPORT_LIMITATIONS = Object.freeze([
    'This report describes deterministic in-silico CDS design checks only.',
    'It does not demonstrate or guarantee expression, yield, folding, biological activity, synthesis acceptance, or regulatory acceptance.',
    'Independent construct review and wet-lab validation are required before experimental reliance.',
]);

const REPORT_CHECK_ACTIONS = Object.freeze({
    cai: 'Compare an alternative implemented method or document why the observed CAI is acceptable for this study.',
    overall_gc: 'Review the host reference band and compare an alternative implemented method before selecting a construct.',
    local_gc: 'Inspect the out-of-range windows and confirm that the local composition is acceptable for the intended construct.',
    type_iis: 'Inspect the listed Type IIS site locations and redesign or explicitly resolve every required site.',
    repeats: 'Inspect repeated regions for synthesis or assembly concerns and document the disposition.',
    homopolymers: 'Inspect the longest homopolymer and confirm it against synthesis requirements.',
    forbidden_motifs: 'Review each detected motif and redesign or document an explicit exception.',
    mfe: 'If RNA structure is decision-relevant, run a fit-for-purpose calculation in an environment with the required dependency.',
});

function finiteNumber(value) {
    const parsed = Number(value);
    return value !== null && value !== '' && Number.isFinite(parsed) ? parsed : null;
}

function reportValue(value, suffix = '') {
    if (value === null || value === undefined || value === '') return 'Not recorded';
    if (typeof value === 'object') return JSON.stringify(value);
    return `${value}${suffix}`;
}

function reportStatusTone(status) {
    if (status === 'FAIL') return 'bad';
    if (status === 'WARNING' || status === 'CONDITIONAL_PASS') return 'warn';
    if (status === 'PASS' || status === 'COMPUTED') return 'good';
    return 'neutral';
}

function checkDetail(criterion, details) {
    if (!details || typeof details !== 'object') return null;
    if (criterion === 'type_iis' && Array.isArray(details.type_iis_sites)) {
        return details.type_iis_sites.map(site => `${site.enzyme || 'site'} at nt ${Number(site.start) + 1}`).join(', ') || null;
    }
    if (criterion === 'forbidden_motifs' && Array.isArray(details.forbidden_motifs)) {
        return details.forbidden_motifs.map(item => typeof item === 'string' ? item : JSON.stringify(item)).join(', ') || null;
    }
    if (criterion === 'local_gc' && details.local_gc_min != null && details.local_gc_max != null) {
        return `Observed range ${details.local_gc_min}–${details.local_gc_max}%`;
    }
    if (criterion === 'repeats' && details.repeat_count != null) return `${details.repeat_count} repeat(s) detected`;
    if (criterion === 'homopolymers' && details.longest_homopolymer != null) return `Longest run: ${details.longest_homopolymer} bp`;
    return null;
}

function normalizedReportChecks(res) {
    const sourceRows = Array.isArray(res.qc_decision_matrix) && res.qc_decision_matrix.length
        ? res.qc_decision_matrix
        : (Array.isArray(res.acceptance_evaluation?.optimized?.criteria)
            ? res.acceptance_evaluation.optimized.criteria
            : []);
    const details = res.acceptance_evaluation?.optimized?.details || {};
    const statusMap = { PASS: 'PASS', FAIL: 'FAIL', WARN: 'WARNING', WARNING: 'WARNING', IGNORED: 'NOT_APPLICABLE' };
    const priority = { FAIL: 0, WARNING: 1, NOT_COMPUTED: 2, NOT_AVAILABLE: 3, NOT_APPLICABLE: 4, PASS: 5, COMPUTED: 5 };
    const checks = sourceRows.map(row => ({
        id: String(row.criterion || 'unknown'),
        label: REPORT_CHECK_LABELS[row.criterion] || String(row.criterion || 'Unknown check'),
        mode: row.mode || 'not_recorded',
        observed: row.observed ?? null,
        threshold: row.threshold ?? null,
        status: statusMap[String(row.result || '').toUpperCase()] || 'NOT_AVAILABLE',
        detail: checkDetail(row.criterion, details),
    }));

    const mfe = getMfeStatus(res);
    checks.push({
        id: 'mfe',
        label: REPORT_CHECK_LABELS.mfe,
        mode: 'informational',
        observed: res.metrics?.mfe_kcal_mol ?? null,
        threshold: 'Informational only',
        status: mfe.status === 'computed' ? 'COMPUTED' : 'NOT_COMPUTED',
        detail: mfe.status === 'computed' ? 'MFE was computed.' : mfe.reason,
    });
    return checks.sort((a, b) => (priority[a.status] ?? 9) - (priority[b.status] ?? 9));
}

function buildResearcherReviewPlan(checks, decision, seed, settingMismatches) {
    const priorities = checks
        .filter(check => ['FAIL', 'WARNING', 'NOT_COMPUTED', 'NOT_AVAILABLE'].includes(check.status))
        .map(check => ({
            check_id: check.id,
            label: check.label,
            status: check.status,
            finding: check.detail || `Observed ${reportValue(check.observed)} against ${reportValue(check.threshold)}.`,
            action: REPORT_CHECK_ACTIONS[check.id] || 'Review this result and document the disposition before downstream use.',
        }));
    settingMismatches.forEach(mismatch => priorities.unshift({
        check_id: 'settings_mismatch',
        label: mismatch.label,
        status: 'FAIL',
        finding: `Requested ${mismatch.requested}; applied ${mismatch.applied}.`,
        action: 'Do not compare or rely on this run until the requested and applied settings are reconciled.',
    }));
    if (seed === null) priorities.push({
        check_id: 'seed',
        label: 'Reproducibility seed',
        status: 'NOT_AVAILABLE',
        finding: 'No seed was recorded for this run.',
        action: 'Set and record a seed before generating a comparison or final candidate.',
    });
    const headline = settingMismatches.length
        ? 'Settings require reconciliation before this run can be interpreted.'
        : decision === 'FAIL'
            ? 'Required checks failed; redesign or document a justified resolution.'
            : decision === 'CONDITIONAL_PASS'
                ? 'No required check failed, but review the warnings before selection.'
                : decision === 'PASS'
                    ? 'Active computational checks passed; continue with independent construct review.'
                    : 'A complete automated decision was not recorded; review the evidence before proceeding.';
    return { headline, priorities };
}

function buildResultsReportModel(res, primary, gcTarget) {
    const checks = normalizedReportChecks(res);
    const summary = res.decision_summary || {};
    const decision = ['PASS', 'CONDITIONAL_PASS', 'FAIL'].includes(res.automated_decision)
        ? res.automated_decision
        : 'NOT_AVAILABLE';
    const rawCandidates = Array.isArray(res.report_candidates)
        ? res.report_candidates
        : (Array.isArray(res.candidates) ? res.candidates : []);
    const custom = res.custom_restriction_sites;
    const requestedNames = Array.isArray(custom?.requested)
        ? custom.requested.map(site => site.name).filter(Boolean)
        : [];
    const selected = requestedNames;
    const effectiveTypeIis = res.acceptance_criteria_snapshot?.type_iis?.enzymes || [];
    const normalizedRequested = [...new Set(selected)].sort();
    const normalizedEffective = [...new Set(effectiveTypeIis)].sort();
    const settingMismatches = normalizedRequested.length && normalizedEffective.length
        && JSON.stringify(normalizedRequested) !== JSON.stringify(normalizedEffective)
        ? [{
            label: 'Type IIS settings mismatch',
            requested: normalizedRequested.join(', '),
            applied: normalizedEffective.join(', '),
        }]
        : [];
    const inputType = res.input_type || res.validation?.input_type || res.provenance?.normalized_input_type || 'not_recorded';
    const outputSequence = primary.optimized_sequence || '';
    const hasComparableInput = inputType === 'cds' && Boolean(state.sequence);
    let baseChanges = null;
    if (hasComparableInput) {
        const comparableLength = Math.min(state.sequence.length, outputSequence.length);
        baseChanges = Math.abs(state.sequence.length - outputSequence.length);
        for (let index = 0; index < comparableLength; index += 1) {
            if (state.sequence[index] !== outputSequence[index]) baseChanges += 1;
        }
    }
    const cai = finiteNumber(primary.metrics?.cai);
    const gc = finiteNumber(primary.metrics?.gc_percent);
    const mfe = getMfeStatus(res);
    const referenceSequence = hasComparableInput ? state.sequence : '';
    const sequenceGcPercent = sequence => {
        if (!sequence) return null;
        const gcCount = (sequence.match(/[GC]/gi) || []).length;
        return (gcCount / sequence.length) * 100;
    };
    const firstWindow = sequence => sequence ? sequence.slice(0, 30) : '';
    const contract = res.design_contract || null;
    const criteriaRows = checks.filter(check => check.id !== 'mfe');
    const requiredFailures = Number.isInteger(summary.required_failure_count)
        ? summary.required_failure_count
        : (criteriaRows.length ? criteriaRows.filter(check => check.status === 'FAIL').length : null);
    const preferredWarnings = Number.isInteger(summary.preferred_warning_count)
        ? summary.preferred_warning_count
        : (criteriaRows.length ? criteriaRows.filter(check => check.status === 'WARNING').length : null);

    const reviewPlan = buildResearcherReviewPlan(checks, decision, res.seed ?? null, settingMismatches);
    return {
        report_schema_version: '1.1',
        identity: {
            result_id: res.result_identifier || res.construct_id || null,
            construct_id: res.construct_id || null,
            result_created_at: res.created_at || null,
            report_generated_at: new Date().toISOString(),
        },
        scientific_review: res.scientific_review || null,
        candidate_scientific_reviews: res.candidate_scientific_reviews || (res.candidates || []).map(c=>({id:c.id,label:c.label,scientific_review:c.scientific_review || null})),
        context: {
            input_type: inputType,
            input_length: finiteNumber(res.original_length ?? res.input_summary?.length),
            host_profile: getResultHostProfile(res),
            profile: res.profile || state.objective || null,
            objective: res.cds_design?.objective || res.profile || state.objective || null,
            engine: res.cds_design?.engine || res.mode || null,
            seed: res.seed ?? null,
        },
        disposition: {
            automated_decision: decision,
            required_failure_count: requiredFailures,
            preferred_warning_count: preferredWarnings,
            unavailable_check_count: checks.filter(check => ['NOT_COMPUTED', 'NOT_AVAILABLE'].includes(check.status)).length,
            explanation: summary.explanation || (decision === 'NOT_AVAILABLE' ? 'Acceptance criteria were not recorded by this result.' : null),
        },
        metrics: {
            cai,
            gc_percent: gc,
            gc_target_min_percent: finiteNumber(gcTarget?.min),
            gc_target_max_percent: finiteNumber(gcTarget?.max),
            mfe_kcal_mol: finiteNumber(res.metrics?.mfe_kcal_mol),
            mfe_status: mfe.status,
            mfe_status_reason: mfe.reason,
            mfe_used: res.metrics?.mfe_used ?? null,
        },
        design_contract: contract ? {
            engine_id: contract.engine_id || res.provenance?.engine_id || null,
            engine_version: contract.engine_version || res.provenance?.engine_version || null,
            engine_status: res.provenance?.engine_status || null,
            scientific_axes: Array.isArray(contract.scientific_axes) ? contract.scientific_axes : [],
            independent_evaluation: contract.independent_evaluation || null,
        } : null,
        comparison: {
            reference_available: hasComparableInput,
            reference_label: hasComparableInput ? 'Input CDS (reference)' : 'Reference CDS not provided',
            candidate_label: contract?.engine_id === 'dp_v2_1_1' ? 'DP v2.1.1 candidate' : 'FactorForge candidate',
            reference_gc_percent: sequenceGcPercent(referenceSequence),
            candidate_gc_percent: gc ?? sequenceGcPercent(outputSequence),
            reference_first_30nt_gc_percent: sequenceGcPercent(firstWindow(referenceSequence)),
            candidate_first_30nt_gc_percent: sequenceGcPercent(firstWindow(outputSequence)),
        },
        checks,
        candidates: rawCandidates.map(candidate => ({
            id: candidate.id || null,
            label: candidate.label || candidate.id || 'Unnamed candidate',
            automated_decision: candidate.automated_decision || 'NOT_AVAILABLE',
            cai: finiteNumber(candidate.cai),
            gc_percent: finiteNumber(candidate.gc_percent),
            required_failure_count: Number.isInteger(candidate.required_failure_count) ? candidate.required_failure_count : null,
            preferred_warning_count: Number.isInteger(candidate.preferred_warning_count) ? candidate.preferred_warning_count : null,
        })),
        sequence_summary: {
            output_length_nt: outputSequence.length || finiteNumber(res.optimized_length),
            amino_acid_identity: finiteNumber(res.constraint_report?.aa_identity),
            nucleotide_changes: baseChanges,
            comparison_available: hasComparableInput,
        },
        process: {
            type_iis_requested: selected,
            type_iis_applied: effectiveTypeIis,
            setting_mismatches: settingMismatches,
            domestication_attempted: Boolean(custom),
            restriction_sites_removed_count: Array.isArray(custom?.removed) ? custom.removed.length : null,
            restriction_sites_unresolved_count: Array.isArray(custom?.unresolved) ? custom.unresolved.length : null,
        },
        provenance: {
            product_version: res.product_version || res.metadata?.product_version || res.cds_design?.product_version || null,
            codon_reference_id: res.codon_reference_id || res.metadata?.codon_reference_id || res.cds_design?.codon_reference_id || null,
            reference_policy_version: res.reference_policy_version || res.metadata?.reference_policy_version || res.cds_design?.reference_policy_version || null,
            gc_reference_band: res.gc_reference_band || res.metadata?.gc_reference_band || res.cds_design?.gc_reference_band || null,
            input_sequence_hash: res.provenance?.input_sequence_hash || null,
            output_cds_hash: res.provenance?.output_cds_hash || null,
            parameter_hash: res.provenance?.parameter_hash || null,
            acceptance_criteria_snapshot: res.acceptance_criteria_snapshot || res.provenance?.acceptance_criteria_snapshot || null,
        },
        interpretation: {
            scope: 'Deterministic in-silico CDS design and pre-synthesis review.',
            headline: reviewPlan.headline,
            review_priorities: reviewPlan.priorities,
            next_steps: ['Resolve and document the review priorities above.', 'Confirm the sequence in its complete construct and assembly context.', 'Perform fit-for-purpose wet-lab validation before experimental reliance.'],
            limitations: [...REPORT_LIMITATIONS],
        },
        structure_viewer: buildStructureViewerPayload(res, primary),
        artifacts: { optimized_sequence: outputSequence, reference_sequence: referenceSequence },
    };
}

function reportPercent(value, digits = 2) {
    return value == null ? 'Not recorded' : `${value.toFixed(digits)}%`;
}

function designContractReportHtml(model) {
    if (!model.design_contract) return '';
    const labels = {
        assembly_feasibility: 'Assembly feasibility',
        codon_adaptation: 'Codon adaptation',
        five_prime_initiation: "5′ initiation-aware scoring",
    };
    const axes = model.design_contract.scientific_axes.map(axis => `<article class="rounded-xl border border-slate-200 bg-white p-3 dark:border-slate-700 dark:bg-slate-950"><p class="font-bold text-slate-900 dark:text-white">${escapeHtml(labels[axis.id] || axis.id)}</p><p class="mt-1 text-[9px] font-extrabold uppercase tracking-widest text-teal-700 dark:text-teal-300">${escapeHtml(axis.evidence_class)}</p></article>`).join('');
    const independent = model.design_contract.independent_evaluation;
    return `<section aria-labelledby="report-contract-title"><div class="flex flex-wrap items-end justify-between gap-2"><div><p class="text-[9px] font-extrabold uppercase tracking-widest text-teal-600 dark:text-teal-300">Declared computational scope</p><h4 id="report-contract-title" class="mt-1 text-sm font-black text-slate-900 dark:text-white">${escapeHtml(`${reportValue(model.design_contract.engine_id)} ${reportValue(model.design_contract.engine_version)}`)}</h4></div><span class="rounded-full bg-amber-100 px-2.5 py-1 text-[9px] font-extrabold uppercase tracking-wide text-amber-800 dark:bg-amber-900/30 dark:text-amber-200">${escapeHtml(reportValue(model.design_contract.engine_status))}</span></div><div class="mt-3 grid gap-2 sm:grid-cols-3">${axes}</div>${independent ? `<div class="mt-2 rounded-xl border border-dashed border-slate-300 p-3 dark:border-slate-600"><p class="font-bold">RNA folding · ${escapeHtml(reportValue(independent.evidence_class))}</p><p class="mt-1 text-[10px] text-slate-500 dark:text-slate-400">Generation status: ${escapeHtml(reportValue(independent.rna_folding))}. This evidence remains separate from candidate generation.</p></div>` : ''}</section>`;
}

function designComparisonReportHtml(model) {
    const reference = model.comparison.reference_available;
    const rows = [
        ['Codon adaptation index (CAI)', 'Not evaluated', model.metrics.cai == null ? 'Not recorded' : model.metrics.cai.toFixed(4), 'Computed metric'],
        ['Global GC', reference ? reportPercent(model.comparison.reference_gc_percent) : 'Not provided', reportPercent(model.comparison.candidate_gc_percent), 'Computed metric'],
        ['First 30 nt GC', reference ? reportPercent(model.comparison.reference_first_30nt_gc_percent) : 'Not provided', reportPercent(model.comparison.candidate_first_30nt_gc_percent), 'Descriptive metric'],
        ['Amino-acid identity', reference ? '100.00%' : 'Input protein', model.sequence_summary.amino_acid_identity == null ? 'Not recorded' : reportPercent(model.sequence_summary.amino_acid_identity * 100), 'Sequence integrity'],
        ['Nucleotide changes', reference ? '0 bp' : 'Not applicable', reference ? `${reportValue(model.sequence_summary.nucleotide_changes)} bp` : 'Not comparable', 'Descriptive metric'],
        ['RNA folding', 'Not evaluated', model.metrics.mfe_status === 'computed' ? reportValue(model.metrics.mfe_kcal_mol, ' kcal/mol') : 'Not computed', 'Independent evaluation'],
    ];
    return `<section aria-labelledby="report-comparison-title"><div class="flex flex-wrap items-end justify-between gap-2"><div><p class="text-[9px] font-extrabold uppercase tracking-widest text-cyan-600 dark:text-cyan-300">Metric-level comparison</p><h4 id="report-comparison-title" class="mt-1 text-sm font-black text-slate-900 dark:text-white">Design Comparison</h4></div><span class="text-[9px] font-bold text-slate-500 dark:text-slate-400">${escapeHtml(reference ? 'Reference CDS available' : 'No reference CDS · candidate-only view')}</span></div><div class="mt-3 overflow-x-auto rounded-xl border border-slate-200 dark:border-slate-700"><table class="w-full min-w-[620px] text-left text-[10px]"><thead class="bg-slate-100 text-slate-600 dark:bg-slate-800 dark:text-slate-300"><tr><th class="p-2.5">Metric</th><th class="p-2.5">${escapeHtml(model.comparison.reference_label)}</th><th class="p-2.5">${escapeHtml(model.comparison.candidate_label)}</th><th class="p-2.5">Evidence</th></tr></thead><tbody class="divide-y divide-slate-200 bg-white dark:divide-slate-700 dark:bg-slate-950">${rows.map(row => `<tr><th scope="row" class="p-2.5 font-semibold">${escapeHtml(row[0])}</th><td class="p-2.5 font-mono">${escapeHtml(row[1])}</td><td class="p-2.5 font-mono font-bold text-cyan-700 dark:text-cyan-300">${escapeHtml(row[2])}</td><td class="p-2.5 text-slate-500 dark:text-slate-400">${escapeHtml(row[3])}</td></tr>`).join('')}</tbody></table></div><p class="mt-2 text-[10px] text-slate-500 dark:text-slate-400">Differences describe deterministic in-silico metrics and do not establish expression, yield, or biological superiority.</p></section>`;
}

function reportCheckRowsHtml(checks) {
    return checks.map(check => `<article class="rounded-xl border border-slate-200 dark:border-slate-700 bg-white dark:bg-slate-950 p-3">
        <div class="flex items-start justify-between gap-3">
            <h5 class="font-bold text-slate-800 dark:text-slate-100">${escapeHtml(check.label)}</h5>
            <span class="shrink-0 text-[9px] font-extrabold tracking-wide" data-report-status="${escapeHtml(check.status)}">${escapeHtml(check.status.replaceAll('_', ' '))}</span>
        </div>
        <dl class="mt-2 grid grid-cols-2 gap-2 text-[10px]">
            <div><dt class="uppercase tracking-widest text-slate-400">Observed</dt><dd class="mt-0.5 font-mono break-all">${escapeHtml(reportValue(check.observed))}</dd></div>
            <div><dt class="uppercase tracking-widest text-slate-400">Policy</dt><dd class="mt-0.5">${escapeHtml(`${check.mode} · ${reportValue(check.threshold)}`)}</dd></div>
        </dl>
        ${check.detail ? `<p class="mt-2 text-[10px] text-slate-500 dark:text-slate-400">${escapeHtml(check.detail)}</p>` : ''}
    </article>`).join('');
}

function reportProvenanceHtml(model) {
    const fields = [
        ['Product version', model.provenance.product_version],
        ['Codon reference', model.provenance.codon_reference_id],
        ['Reference policy', model.provenance.reference_policy_version],
        ['GC reference band', model.provenance.gc_reference_band],
        ['Seed', model.context.seed === null ? 'Not specified' : model.context.seed],
        ['Result created (UTC)', model.identity.result_created_at],
        ['Input SHA-256', model.provenance.input_sequence_hash],
        ['Output SHA-256', model.provenance.output_cds_hash],
        ['Parameter SHA-256', model.provenance.parameter_hash],
    ];
    return fields.map(([label, value]) => `<div class="min-w-0"><dt class="text-[9px] uppercase tracking-widest text-slate-500 dark:text-slate-400">${escapeHtml(label)}</dt><dd class="mt-0.5 font-mono break-all text-[10px]">${escapeHtml(reportValue(value))}</dd></div>`).join('');
}

function candidateReportHtml(candidates) {
    if (candidates.length <= 1) return '';
    return `<section aria-labelledby="report-candidates-title" class="pt-2">
        <h4 id="report-candidates-title" class="text-[10px] font-extrabold uppercase tracking-widest text-slate-500 dark:text-slate-400 mb-2">Candidate comparison</h4>
        <div class="overflow-x-auto"><table class="w-full text-left text-xs">
            <thead><tr class="text-slate-500 dark:text-slate-400"><th class="py-1 pr-3">Candidate</th><th class="py-1 pr-3">Decision</th><th class="py-1 pr-3">CAI</th><th class="py-1 pr-3">GC%</th><th class="py-1">Issues</th></tr></thead>
            <tbody class="divide-y divide-slate-200 dark:divide-slate-700">${candidates.map(candidate => `<tr><th scope="row" class="py-2 pr-3 font-semibold">${escapeHtml(candidate.label)}</th><td class="py-2 pr-3 font-bold">${escapeHtml(candidate.automated_decision.replaceAll('_', ' '))}</td><td class="py-2 pr-3 font-mono">${escapeHtml(candidate.cai == null ? 'Not recorded' : candidate.cai.toFixed(3))}</td><td class="py-2 pr-3 font-mono">${escapeHtml(candidate.gc_percent == null ? 'Not recorded' : candidate.gc_percent.toFixed(1))}</td><td class="py-2">${escapeHtml(`${reportValue(candidate.required_failure_count)} fail · ${reportValue(candidate.preferred_warning_count)} warn`)}</td></tr>`).join('')}</tbody>
        </table></div>
    </section>`;
}

function reportPriorityRowsHtml(priorities) {
    if (!priorities.length) return '<p class="text-xs text-slate-500 dark:text-slate-400">No unresolved computational review priority was recorded.</p>';
    return priorities.map((item, index) => `<article class="rounded-xl border border-slate-200 dark:border-slate-700 bg-white dark:bg-slate-950 p-3" data-review-priority="${escapeHtml(item.check_id)}">
        <div class="flex items-start gap-3"><span class="flex h-6 w-6 shrink-0 items-center justify-center rounded-full bg-slate-100 dark:bg-slate-800 text-[10px] font-black">${index + 1}</span><div class="min-w-0 flex-1">
            <div class="flex items-start justify-between gap-2"><h5 class="font-bold text-slate-900 dark:text-white">${escapeHtml(item.label)}</h5><span class="text-[9px] font-extrabold" data-report-status="${escapeHtml(item.status)}">${escapeHtml(item.status.replaceAll('_', ' '))}</span></div>
            <p class="mt-1 text-[10px] text-slate-500 dark:text-slate-400">${escapeHtml(item.finding)}</p>
            <p class="mt-2 text-[10px] font-semibold text-slate-700 dark:text-slate-200"><span class="uppercase tracking-wider text-slate-400">Next:</span> ${escapeHtml(item.action)}</p>
        </div></div>
    </article>`).join('');
}

function reportSettingsHtml(model) {
    const requested = model.process.type_iis_requested.length ? model.process.type_iis_requested.join(', ') : 'None recorded';
    const applied = model.process.type_iis_applied.length ? model.process.type_iis_applied.join(', ') : 'Not recorded';
    const mismatch = model.process.setting_mismatches.length > 0;
    const fields = [
        ['Host', model.context.host_profile, model.context.host_profile, false],
        ['Method / profile', model.context.objective, model.context.profile, false],
        ['Type IIS enzymes', requested, applied, mismatch],
        ['Seed', model.context.seed === null ? 'Not specified' : model.context.seed, model.context.seed === null ? 'Not recorded' : model.context.seed, false],
        ['GC reference', model.provenance.gc_reference_band, model.provenance.gc_reference_band, false],
    ];
    return fields.map(([label, requestedValue, appliedValue, differs]) => `<tr class="border-b border-slate-200 dark:border-slate-700 ${differs ? 'bg-rose-50 dark:bg-rose-900/20' : ''}"><th scope="row" class="py-2 pr-2 text-left font-semibold">${escapeHtml(label)}</th><td class="py-2 pr-2">${escapeHtml(reportValue(requestedValue))}</td><td class="py-2">${escapeHtml(reportValue(appliedValue))}${differs ? ' <b class="text-rose-700 dark:text-rose-300">Mismatch</b>' : ''}</td></tr>`).join('');
}

function renderResultsReport(res, primary, gcTarget, reportData = null) {
    if (!elements.resultsReport || !elements.resultsReportBody) return;
    const model = reportData || buildResultsReportModel(res, primary, gcTarget);
    const disposition = model.disposition.automated_decision;
    const comparisonText = model.context.input_type === 'protein'
        ? `Protein input · amino-acid identity ${model.sequence_summary.amino_acid_identity == null ? 'Not recorded' : `${(model.sequence_summary.amino_acid_identity * 100).toFixed(2)}%`}`
        : `CDS input · nucleotide changes ${reportValue(model.sequence_summary.nucleotide_changes)}`;
    const domesticationText = model.process.domestication_attempted
        ? `Attempted · ${reportValue(model.process.restriction_sites_removed_count)} removed · ${reportValue(model.process.restriction_sites_unresolved_count)} unresolved`
        : 'Not attempted';

    elements.resultsReportBody.innerHTML = `<article aria-labelledby="design-review-report-title" class="space-y-4">
        <header>
            <h3 id="design-review-report-title" class="text-sm font-black text-slate-900 dark:text-white">Researcher Decision Report</h3>
            <p class="mt-1 text-[10px] text-slate-500 dark:text-slate-400">${escapeHtml(`${reportValue(model.identity.result_id)} · ${model.context.host_profile} · ${model.context.profile} · ${model.context.seed === null ? 'seed not specified' : `seed=${model.context.seed}`}`)}</p>
        </header>
        <section aria-labelledby="report-outcome-title" class="rounded-2xl border p-4 ${disposition === 'FAIL' ? 'border-rose-300 bg-rose-50 dark:border-rose-800 dark:bg-rose-900/20' : disposition === 'CONDITIONAL_PASS' ? 'border-amber-300 bg-amber-50 dark:border-amber-800 dark:bg-amber-900/20' : 'border-emerald-300 bg-emerald-50 dark:border-emerald-800 dark:bg-emerald-900/20'}">
            <p class="text-[9px] font-extrabold uppercase tracking-widest opacity-70">Researcher decision brief</p>
            <div class="mt-2 flex flex-wrap items-end justify-between gap-2"><h4 id="report-outcome-title" class="text-xl font-black">${escapeHtml(disposition.replaceAll('_', ' '))}</h4><p class="text-[10px] font-bold">${escapeHtml(`${reportValue(model.disposition.required_failure_count)} required fail · ${reportValue(model.disposition.preferred_warning_count)} warning · ${reportValue(model.disposition.unavailable_check_count)} unavailable`)}</p></div>
            <p class="mt-2 text-xs font-semibold">${escapeHtml(model.interpretation.headline)}</p>
            <p class="mt-1 text-[10px] opacity-80">${escapeHtml(reportValue(model.disposition.explanation))}</p>
        </section>
        ${designContractReportHtml(model)}
        ${scientificReviewHtml(model)}${designComparisonReportHtml(model)}
        <section aria-labelledby="report-priorities-title">
            <h4 id="report-priorities-title" class="text-[10px] font-extrabold uppercase tracking-widest text-slate-500 dark:text-slate-400 mb-2">Review priorities and next actions</h4>
            <div class="grid gap-2">${reportPriorityRowsHtml(model.interpretation.review_priorities)}</div>
        </section>
        <section aria-labelledby="report-settings-title"><h4 id="report-settings-title" class="text-[10px] font-extrabold uppercase tracking-widest text-slate-500 dark:text-slate-400 mb-2">Requested vs applied settings</h4><div class="overflow-x-auto"><table class="w-full text-[10px]"><thead><tr class="text-left text-slate-500"><th class="pb-2">Setting</th><th class="pb-2">Requested</th><th class="pb-2">Applied / recorded</th></tr></thead><tbody>${reportSettingsHtml(model)}</tbody></table></div></section>
        <section aria-labelledby="report-sequence-title" class="rounded-xl border border-slate-200 dark:border-slate-700 p-3">
            <h4 id="report-sequence-title" class="text-[10px] font-extrabold uppercase tracking-widest text-slate-500 dark:text-slate-400">Sequence and process</h4>
            <p class="mt-2 font-semibold">${escapeHtml(comparisonText)}</p>
            <p class="mt-1">Output length: ${escapeHtml(reportValue(model.sequence_summary.output_length_nt, ' nt'))}</p>
            <p class="mt-1">Type IIS requested: ${escapeHtml(model.process.type_iis_requested.length ? model.process.type_iis_requested.join(', ') : 'None recorded')}</p>
            <p class="mt-1">Domestication: ${escapeHtml(domesticationText)}</p>
            <p class="mt-1">MFE: ${escapeHtml(model.metrics.mfe_status === 'computed' ? `Computed${model.metrics.mfe_kcal_mol == null ? '' : ` · ${model.metrics.mfe_kcal_mol} kcal/mol`}` : `Not computed · ${model.metrics.mfe_status_reason}`)}</p>
        </section>
        ${candidateReportHtml(model.candidates)}
        <details class="rounded-xl border border-slate-200 dark:border-slate-700 p-3"><summary class="cursor-pointer text-[10px] font-extrabold uppercase tracking-widest">All computational checks (${model.checks.length})</summary><div class="mt-3 grid gap-2">${reportCheckRowsHtml(model.checks)}</div></details>
        <details class="rounded-xl border border-slate-200 dark:border-slate-700 p-3"><summary class="cursor-pointer text-[10px] font-extrabold uppercase tracking-widest">Reproducibility and provenance</summary><dl class="mt-3 grid grid-cols-1 sm:grid-cols-2 gap-3">${reportProvenanceHtml(model)}</dl></details>
        <section aria-labelledby="report-interpretation-title" class="grid gap-2">
            <div class="rounded-xl bg-blue-50 dark:bg-blue-900/20 p-3"><h4 id="report-interpretation-title" class="font-extrabold">Interpretation</h4><p class="mt-1">${escapeHtml(model.interpretation.scope)}</p></div>
            <div class="rounded-xl bg-emerald-50 dark:bg-emerald-900/20 p-3"><h4 class="font-extrabold">Downstream handoff checklist</h4><ul class="mt-1 list-disc pl-4">${model.interpretation.next_steps.map(item => `<li>${escapeHtml(item)}</li>`).join('')}</ul></div>
            <div class="rounded-xl bg-amber-50 dark:bg-amber-900/20 p-3"><h4 class="font-extrabold">Limitations</h4><ul class="mt-1 list-disc pl-4">${model.interpretation.limitations.map(item => `<li>${escapeHtml(item)}</li>`).join('')}</ul></div>
        </section>
        <p class="rounded-xl border border-amber-200 dark:border-amber-800 p-3 text-[10px] text-amber-800 dark:text-amber-200">The HTML report contains the optimized DNA sequence. Treat it according to your sequence-data policy. The evidence JSON intentionally excludes raw sequences.</p>
        <div class="grid grid-cols-1 sm:grid-cols-2 gap-2">
            <button type="button" id="downloadResultsReportBtn" class="w-full rounded-xl bg-slate-800 hover:bg-slate-700 dark:bg-slate-700 dark:hover:bg-slate-600 text-white text-xs font-bold py-2.5 transition-colors">📄 Download Report (HTML)</button>
            <button type="button" id="downloadEvidenceRecordBtn" class="w-full rounded-xl bg-indigo-700 hover:bg-indigo-600 text-white text-xs font-bold py-2.5 transition-colors">🧾 Download Evidence Record (JSON)</button>
        </div>
    </article>`;
    elements.resultsReport.classList.remove('hidden');
    document.getElementById('downloadResultsReportBtn')?.addEventListener('click', () => downloadResultsReportHtml(model));
    document.getElementById('downloadEvidenceRecordBtn')?.addEventListener('click', () => downloadEvidenceRecordJson(model));
}

function reportFileStem(model) {
    const raw = model.identity.result_id || model.identity.report_generated_at || Date.now();
    return String(raw).replace(/[^A-Za-z0-9._-]+/g, '-').replace(/^-+|-+$/g, '') || String(Date.now());
}

function downloadTextArtifact(content, mimeType, fileName) {
    const blob = new Blob([content], { type: mimeType });
    const url = window.URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = fileName;
    link.click();
    window.URL.revokeObjectURL(url);
}

function evidenceRecordFromModel(model) {
    const { artifacts: _artifacts, structure_viewer: _structureViewer, ...evidence } = model;
    return evidence;
}

function standaloneStructure3dReportSection(model) {
    const sv = model.structure_viewer;
    if (!sv || !sv.cartoonPdb || !Array.isArray(sv.residues) || !sv.residues.length) return { html: '', script: '' };
    const is3dp = sv.mode === '3dp';
    const safeSv = JSON.stringify(sv).replace(/</g, '\\u003c');
    const titleText = is3dp
        ? 'Original Structure Coordinates vs 3DP Structure-Aware Codon Mapping'
        : 'Original Structure Coordinates vs 1DP Baseline Codon Mapping';
    const subText = is3dp
        ? `Structure source: ${escapeHtml(sv.source_label)} · Left: input B-factor/confidence values · Right: 3DP translational pause targets (${sv.pauseCount} residues)`
        : `Structure source: ${escapeHtml(sv.source_label)} · Left: input B-factor/confidence values · Right: 1DP preferred codon substitutions`;
    const btns = is3dp ? `
        <button type="button" data-r3d-pair="orig_3dp" style="padding:5px 10px;border-radius:8px;border:none;background:#059669;color:#fff;font-size:11px;font-weight:800;cursor:pointer;">Original | 3DP</button>
        <button type="button" data-r3d-pair="orig_1dp" style="padding:5px 10px;border-radius:8px;border:none;background:#1e293b;color:#cbd5e1;font-size:11px;font-weight:700;cursor:pointer;">Original | 1DP</button>
        <button type="button" data-r3d-pair="1dp_3dp" style="padding:5px 10px;border-radius:8px;border:none;background:#1e293b;color:#cbd5e1;font-size:11px;font-weight:700;cursor:pointer;">1DP | 3DP</button>
    ` : `<span style="padding:4px 9px;border-radius:8px;background:#0c4a6e;color:#bae6fd;font-size:10px;font-weight:800;">Original | 1DP Side-by-Side</span>`;
    const html = `<section><h2>3D Codon Structure Viewer (${is3dp ? 'Original | 3DP' : 'Original | 1DP'})</h2><div style="background:#0f172a;color:#f8fafc;border-radius:16px;padding:16px;border:1px solid #1e293b;"><div style="display:flex;flex-wrap:wrap;justify-content:space-between;align-items:center;gap:10px;margin-bottom:12px;"><div><strong style="font-size:15px;color:#fff;">${titleText}</strong><div style="font-size:11px;color:#94a3b8;margin-top:2px;">${subText}</div></div><div style="display:flex;flex-wrap:wrap;gap:6px;align-items:center;">${btns}<button type="button" id="r3dReset" style="padding:5px 10px;border-radius:8px;border:1px solid #334155;background:#1e293b;color:#e2e8f0;font-size:11px;font-weight:700;cursor:pointer;">Reset View</button><button type="button" id="r3dSpin" style="padding:5px 10px;border-radius:8px;border:1px solid #334155;background:#1e293b;color:#e2e8f0;font-size:11px;font-weight:700;cursor:pointer;">Auto Spin: Off</button></div></div><div id="r3dGrid" style="display:grid;grid-template-columns:1fr 1fr;gap:10px;"><div id="r3dPaneLeft" style="position:relative;height:380px;border-radius:12px;overflow:hidden;border:1px solid #1e293b;background:#070b14;"><div style="position:absolute;top:10px;left:10px;z-index:2;background:rgba(15,23,42,0.88);border:1px solid #334155;border-radius:8px;padding:6px 10px;pointer-events:none;"><b id="r3dLeftTitle" style="color:#38bdf8;font-size:12px;">LEFT · Original Structure (Input Structure)</b><div id="r3dLeftSub" style="font-size:10px;color:#94a3b8;">Input coordinates (B-factor/confidence) · White sphere = pause target</div></div><canvas id="r3dCanvasLeft" style="width:100%;height:100%;display:block;cursor:grab;"></canvas></div><div id="r3dPaneRight" style="position:relative;height:380px;border-radius:12px;overflow:hidden;border:1px solid #1e293b;background:#070b14;"><div style="position:absolute;top:10px;left:10px;z-index:2;background:rgba(15,23,42,0.88);border:1px solid #334155;border-radius:8px;padding:6px 10px;pointer-events:none;"><b id="r3dRightTitle" style="color:#34d399;font-size:12px;">${is3dp ? 'RIGHT · 3DP Structure-Aware Design (Post-CDS)' : 'RIGHT · 1DP Baseline Design (Post-CDS)'}</b><div id="r3dRightSub" style="font-size:10px;color:#94a3b8;">CAI ${(is3dp ? (sv.cai3dp ?? sv.cai1dp) : sv.cai1dp).toFixed(3)} · GC ${(is3dp ? (sv.gc3dp ?? sv.gc1dp) : sv.gc1dp).toFixed(1)}%</div></div><canvas id="r3dCanvasRight" style="width:100%;height:100%;display:block;cursor:grab;"></canvas></div></div><div style="margin-top:10px;display:flex;flex-wrap:wrap;justify-content:space-between;align-items:center;gap:8px;background:#070b14;border:1px solid #1e293b;border-radius:10px;padding:8px 12px;font-size:11px;color:#cbd5e1;"><div style="display:flex;flex-wrap:wrap;gap:12px;"><span>🔵 Left: Input B-factor/confidence</span><span>🔴 Right Red: Pause-site slowdown</span><span>🔵 Right Cyan: Synonymous shift</span><span>⚫ Right Slate: Unchanged codon</span></div><div id="r3dHover" style="font-family:monospace;color:#6ee7b7;">Click or hover residue for codon details</div></div></div></section>`;
    const script = `<script>(()=>{const sv=${safeSv},rs=sv.residues||[];if(!rs.length)return;const cx=rs.reduce((s,r)=>s+r.x,0)/rs.length,cy=rs.reduce((s,r)=>s+r.y,0)/rs.length,cz=rs.reduce((s,r)=>s+r.z,0)/rs.length,pts=rs.map(r=>({...r,cx:r.x-cx,cy:r.y-cy,cz:r.z-cz}));let lm='orig',rm=sv.mode==='3dp'?'3dp':'1dp',sel=sv.defaultSelectedPos||1,cmL={rx:0.22,ry:-0.55,zoom:1},cmR={rx:0.22,ry:-0.55,zoom:1},spin=false;const hEl=document.getElementById('r3dHover');const cr=(p0,p1,p2,p3,t)=>{const t2=t*t,t3=t2*t;return[0.5*((2*p1[0])+(-p0[0]+p2[0])*t+(2*p0[0]-5*p1[0]+4*p2[0]-p3[0])*t2+(-p0[0]+3*p1[0]-3*p2[0]+p3[0])*t3),0.5*((2*p1[1])+(-p0[1]+p2[1])*t+(2*p0[1]-5*p1[1]+4*p2[1]-p3[1])*t2+(-p0[1]+3*p1[1]-3*p2[1]+p3[1])*t3),0.5*((2*p1[2])+(-p0[2]+p2[2])*t+(2*p0[2]-5*p1[2]+4*p2[2]-p3[2])*t2+(-p0[2]+3*p1[2]-3*p2[2]+p3[2])*t3)]};const sh=(hex,f)=>{const h=String(hex||'#38bdf8').replace('#','');return'rgb('+Math.min(255,Math.max(0,Math.round(parseInt(h.slice(0,2),16)*f)))+','+Math.min(255,Math.max(0,Math.round(parseInt(h.slice(2,4),16)*f)))+','+Math.min(255,Math.max(0,Math.round(parseInt(h.slice(4,6),16)*f)))+')'};const cLeft=r=>r.plddt>=90?'#1d4ed8':r.plddt>=70?'#38bdf8':r.plddt>=50?'#eab308':'#f97316';const cRight=(r,m)=>m==='3dp'?((r.isPause&&(r.isDiff||r.isDiffOrig3dp))?'#ef4444':r.isPause?'#eab308':(r.isDiff||r.isDiffOrig3dp)?'#38bdf8':'#334155'):(r.isDiffOrig1dp?(r.freq1dp>=0.35?'#10b981':'#38bdf8'):'#334155');function bind(cv,gM,gC,sC){if(!cv)return()=>{};let dr=false,lx=0,ly=0,pr=[];function proj(pt,w,h,cm){const x1=pt[0]*Math.cos(cm.ry)+pt[2]*Math.sin(cm.ry),z1=-pt[0]*Math.sin(cm.ry)+pt[2]*Math.cos(cm.ry),y2=pt[1]*Math.cos(cm.rx)-z1*Math.sin(cm.rx),z2=pt[1]*Math.sin(cm.rx)+z1*Math.cos(cm.rx),sc=(560*cm.zoom)/(160+z2);return{sx:w*0.5+x1*sc,sy:h*0.5-y2*sc,z:z2,sc}}function draw(){const rc=cv.getBoundingClientRect();if(!rc.width||!rc.height)return;const d=devicePixelRatio||1;cv.width=rc.width*d;cv.height=rc.height*d;const x=cv.getContext('2d');x.setTransform(d,0,0,d,0,0);x.fillStyle='#070b14';x.fillRect(0,0,rc.width,rc.height);const mk=gM(),cm=gC(),segs=[];for(let i=0;i<pts.length-1;i++){const p0=pts[Math.max(0,i-1)],p1=pts[i],p2=pts[i+1],p3=pts[Math.min(pts.length-1,i+2)],c0=[p0.cx,p0.cy,p0.cz],c1=[p1.cx,p1.cy,p1.cz],c2=[p2.cx,p2.cy,p2.cz],c3=[p3.cx,p3.cy,p3.cz],col=mk==='orig'?cLeft(p1):cRight(p1,mk);for(let s=0;s<4;s++){const a=proj(cr(c0,c1,c2,c3,s/4),rc.width,rc.height,cm),b=proj(cr(c0,c1,c2,c3,(s+1)/4),rc.width,rc.height,cm);segs.push({z:(a.z+b.z)*0.5,a,b,col,w:2.0*((a.sc+b.sc)*0.5)})}}pr=pts.map(r=>({...proj([r.cx,r.cy,r.cz],rc.width,rc.height,cm),r}));segs.sort((a,b)=>b.z-a.z).forEach(s=>{const df=Math.max(0.48,Math.min(1.22,1-s.z/95));x.beginPath();x.moveTo(s.a.sx,s.a.sy);x.lineTo(s.b.sx,s.b.sy);x.strokeStyle=sh(s.col,df);x.lineWidth=Math.max(1.8,s.w);x.lineCap='round';x.stroke()});const sp=pr.find(p=>p.r.pos===sel);if(sp){const r=sp.r,tg=mk==='orig'?r.aa+r.pos+' · Orig:'+r.codonOrig+' (pLDDT='+r.plddt.toFixed(1)+')':r.aa+r.pos+' · '+(mk==='3dp'?'3DP:'+r.codon3dp:'1DP:'+r.codon1dp);x.font='bold 11px monospace';const tw=x.measureText(tg).width+12,tx=Math.min(rc.width-tw-6,Math.max(6,sp.sx+8)),ty=Math.max(24,Math.min(rc.height-28,sp.sy-8));x.fillStyle='rgba(15,23,42,0.92)';x.fillRect(tx,ty-14,tw,20);x.fillStyle='#fff';x.fillText(tg,tx+6,ty)}}cv.addEventListener('mousedown',e=>{dr=true;lx=e.clientX;ly=e.clientY});addEventListener('mouseup',()=>{dr=false});addEventListener('mousemove',e=>{if(!dr)return;const cm=gC();sC({rx:Math.max(-1.45,Math.min(1.45,cm.rx+(e.clientY-ly)*0.008)),ry:cm.ry+(e.clientX-lx)*0.008,zoom:cm.zoom});lx=e.clientX;ly=e.clientY;draw()});cv.addEventListener('click',e=>{const rc=cv.getBoundingClientRect(),mx=e.clientX-rc.left,my=e.clientY-rc.top;let b=null,md=18;pr.forEach(p=>{const d=Math.hypot(p.sx-mx,p.sy-my);if(d<md){md=d;b=p}});if(b){sel=b.r.pos;const r=b.r;if(hEl)hEl.textContent='Pos '+r.pos+' ('+r.aa+') | Orig: '+r.codonOrig+' → 1DP: '+r.codon1dp+(r.codon3dp?' → 3DP: '+r.codon3dp:'')+' | pLDDT '+r.plddt.toFixed(1);dL();dR()}});cv.addEventListener('wheel',e=>{e.preventDefault();const cm=gC();sC({...cm,zoom:Math.max(0.45,Math.min(3.5,cm.zoom*(e.deltaY>0?0.92:1.08)))});draw()},{passive:false});return draw}const dL=bind(document.getElementById('r3dCanvasLeft'),()=>lm,()=>cmL,n=>{cmL=n}),dR=bind(document.getElementById('r3dCanvasRight'),()=>rm,()=>cmR,n=>{cmR=n});document.querySelectorAll('[data-r3d-pair]').forEach(btn=>btn.addEventListener('click',()=>{const p=btn.getAttribute('data-r3d-pair');document.querySelectorAll('[data-r3d-pair]').forEach(b=>{b.style.background=b===btn?'#059669':'#1e293b';b.style.color=b===btn?'#fff':'#cbd5e1'});if(p==='orig_3dp'){lm='orig';rm='3dp'}else if(p==='orig_1dp'){lm='orig';rm='1dp'}else{lm='1dp';rm='3dp'}dL();dR()}));document.getElementById('r3dReset')?.addEventListener('click',()=>{cmL={rx:0.22,ry:-0.55,zoom:1};cmR={rx:0.22,ry:-0.55,zoom:1};dL();dR()});const sb=document.getElementById('r3dSpin');sb?.addEventListener('click',()=>{spin=!spin;sb.textContent='Auto Spin: '+(spin?'On':'Off')});function loop(){if(spin){cmL.ry+=0.0045;cmR.ry+=0.0045;dL();dR()}requestAnimationFrame(loop)}addEventListener('resize',()=>{dL();dR()});dL();dR();requestAnimationFrame(loop)})();<\/script>`;
    return { html, script };
}

function standaloneReportHtml(model) {
    const reportTheme = document.documentElement.classList.contains('dark') ? 'dark' : 'light';
    const statusColor = status => ({ good: '#059669', warn: '#d97706', bad: '#e11d48', neutral: '#475569' })[reportStatusTone(status)];
    const decisionColor = statusColor(model.disposition.automated_decision);
    const priorities = model.interpretation.review_priorities.map((item, index) => `<article class="priority" style="border-left-color:${statusColor(item.status)}"><div class="priority-number">${index + 1}</div><div><h3>${escapeHtml(item.label)} <small style="color:${statusColor(item.status)}">${escapeHtml(item.status.replaceAll('_', ' '))}</small></h3><p>${escapeHtml(item.finding)}</p><p><b>Next:</b> ${escapeHtml(item.action)}</p></div></article>`).join('') || '<p>No unresolved computational review priority was recorded.</p>';
    const requestedTypeIis = model.process.type_iis_requested.length ? model.process.type_iis_requested.join(', ') : 'None recorded';
    const appliedTypeIis = model.process.type_iis_applied.length ? model.process.type_iis_applied.join(', ') : 'Not recorded';
    const settingsMismatch = model.process.setting_mismatches.length > 0;
    const settings = [
        ['Host', model.context.host_profile, model.context.host_profile, false],
        ['Method / profile', model.context.objective, model.context.profile, false],
        ['Type IIS enzymes', requestedTypeIis, appliedTypeIis, settingsMismatch],
        ['Seed', model.context.seed === null ? 'Not specified' : model.context.seed, model.context.seed === null ? 'Not recorded' : model.context.seed, false],
        ['GC reference', model.provenance.gc_reference_band, model.provenance.gc_reference_band, false],
    ].map(([label, requested, applied, mismatch]) => `<tr${mismatch ? ' class="mismatch"' : ''}><th scope="row">${escapeHtml(label)}</th><td>${escapeHtml(reportValue(requested))}</td><td>${escapeHtml(reportValue(applied))}${mismatch ? ' <b>Mismatch</b>' : ''}</td></tr>`).join('');
    const checks = model.checks.map(check => `<tr><th scope="row">${escapeHtml(check.label)}</th><td>${escapeHtml(reportValue(check.observed))}</td><td>${escapeHtml(`${check.mode} · ${reportValue(check.threshold)}`)}</td><td><b style="color:${statusColor(check.status)}">${escapeHtml(check.status.replaceAll('_', ' '))}</b>${check.detail ? `<br><small>${escapeHtml(check.detail)}</small>` : ''}</td></tr>`).join('');
    const candidates = model.candidates.length > 1 ? `<section><h2>Candidate comparison</h2><div class="scroll"><table><thead><tr><th>Candidate</th><th>Decision</th><th>CAI</th><th>GC%</th><th>Required failures</th><th>Warnings</th></tr></thead><tbody>${model.candidates.map(candidate => `<tr><th scope="row">${escapeHtml(candidate.label)}</th><td>${escapeHtml(candidate.automated_decision)}</td><td>${escapeHtml(candidate.cai == null ? 'Not recorded' : candidate.cai.toFixed(3))}</td><td>${escapeHtml(candidate.gc_percent == null ? 'Not recorded' : candidate.gc_percent.toFixed(1))}</td><td>${escapeHtml(reportValue(candidate.required_failure_count))}</td><td>${escapeHtml(reportValue(candidate.preferred_warning_count))}</td></tr>`).join('')}</tbody></table></div></section>` : '';
    const axisLabels = { assembly_feasibility: 'Assembly feasibility', codon_adaptation: 'Codon adaptation', five_prime_initiation: "5′ initiation-aware scoring" };
    const standaloneContract = model.design_contract ? `<section><h2>Declared computational scope</h2><div class="section-card"><div class="section-head"><strong>${escapeHtml(`${reportValue(model.design_contract.engine_id)} ${reportValue(model.design_contract.engine_version)}`)}</strong><span class="badge">${escapeHtml(reportValue(model.design_contract.engine_status))}</span></div><div class="axis-grid">${model.design_contract.scientific_axes.map(axis => `<article><b>${escapeHtml(axisLabels[axis.id] || axis.id)}</b><small>${escapeHtml(axis.evidence_class)}</small></article>`).join('')}</div>${model.design_contract.independent_evaluation ? `<p class="independent"><b>RNA folding · ${escapeHtml(reportValue(model.design_contract.independent_evaluation.evidence_class))}</b><br><small>Generation status: ${escapeHtml(reportValue(model.design_contract.independent_evaluation.rna_folding))}. Kept separate from candidate generation.</small></p>` : ''}</div></section>` : '';
    const comparisonRows = [
        ['Codon adaptation index (CAI)', 'Not evaluated', model.metrics.cai == null ? 'Not recorded' : model.metrics.cai.toFixed(4), 'Computed metric'],
        ['Global GC', model.comparison.reference_available ? reportPercent(model.comparison.reference_gc_percent) : 'Not provided', reportPercent(model.comparison.candidate_gc_percent), 'Computed metric'],
        ['First 30 nt GC', model.comparison.reference_available ? reportPercent(model.comparison.reference_first_30nt_gc_percent) : 'Not provided', reportPercent(model.comparison.candidate_first_30nt_gc_percent), 'Descriptive metric'],
        ['Amino-acid identity', model.comparison.reference_available ? '100.00%' : 'Input protein', model.sequence_summary.amino_acid_identity == null ? 'Not recorded' : reportPercent(model.sequence_summary.amino_acid_identity * 100), 'Sequence integrity'],
        ['Nucleotide changes', model.comparison.reference_available ? '0 bp' : 'Not applicable', model.comparison.reference_available ? `${reportValue(model.sequence_summary.nucleotide_changes)} bp` : 'Not comparable', 'Descriptive metric'],
        ['RNA folding', 'Not evaluated', model.metrics.mfe_status === 'computed' ? reportValue(model.metrics.mfe_kcal_mol, ' kcal/mol') : 'Not computed', 'Independent evaluation'],
    ];
    const standaloneComparison = `<section><h2>Design comparison</h2><p class="muted">${escapeHtml(model.comparison.reference_available ? 'Reference CDS available' : 'No reference CDS · candidate-only view')}</p><div class="scroll"><table><thead><tr><th>Metric</th><th>${escapeHtml(model.comparison.reference_label)}</th><th>${escapeHtml(model.comparison.candidate_label)}</th><th>Evidence</th></tr></thead><tbody>${comparisonRows.map(row => `<tr><th scope="row">${escapeHtml(row[0])}</th><td>${escapeHtml(row[1])}</td><td class="candidate-value">${escapeHtml(row[2])}</td><td>${escapeHtml(row[3])}</td></tr>`).join('')}</tbody></table></div><p class="boundary">Differences describe deterministic in-silico metrics and do not establish expression, yield, or biological superiority.</p></section>`;
    const provenanceRows = [
        ['Product version', model.provenance.product_version], ['Engine', model.context.engine], ['Objective', model.context.objective],
        ['Host / profile', `${model.context.host_profile} / ${model.context.profile}`], ['Codon reference', model.provenance.codon_reference_id],
        ['Reference policy', model.provenance.reference_policy_version], ['GC reference band', model.provenance.gc_reference_band],
        ['Seed', model.context.seed === null ? 'Not specified' : model.context.seed], ['Result created (UTC)', model.identity.result_created_at],
        ['Input SHA-256', model.provenance.input_sequence_hash], ['Output SHA-256', model.provenance.output_cds_hash], ['Parameter SHA-256', model.provenance.parameter_hash],
    ].map(([label, value]) => `<dt>${escapeHtml(label)}</dt><dd>${escapeHtml(reportValue(value))}</dd>`).join('');
    const sequence = model.artifacts.optimized_sequence.match(/.{1,60}/g)?.join('\n') || model.artifacts.optimized_sequence;
    const referenceSequence = model.artifacts.reference_sequence?.match(/.{1,60}/g)?.join('\n') || model.artifacts.reference_sequence;
    const structure3dReport = standaloneStructure3dReportSection(model);
    return `<!doctype html><html lang="en" data-theme="${reportTheme}"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>FactorForge Researcher Decision Report — ${escapeHtml(reportValue(model.identity.result_id))}</title><style>
:root{color-scheme:light}body{margin:0;background:#f1f5f9;color:#0f172a;font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}.page{max-width:920px;margin:auto;padding:42px 24px 72px}h1{font-size:30px;margin:4px 0}h2{font-size:14px;text-transform:uppercase;letter-spacing:.08em;color:#475569;margin:28px 0 10px}h3{margin:0}.eyebrow{color:#047857;font-weight:800;text-transform:uppercase;letter-spacing:.08em}.muted,small{color:#64748b}.outcome,.section-card{background:#fff;border:2px solid;border-radius:16px;padding:20px}.outcome strong{display:block;font-size:28px}.counts{font-weight:700}.callout{padding:14px;border-radius:12px;background:#fff7ed;border:1px solid #fdba74}.section-head{display:flex;justify-content:space-between;gap:12px;align-items:center}.badge{border-radius:999px;background:#fef3c7;color:#92400e;padding:5px 9px;font-size:10px;font-weight:800;text-transform:uppercase}.axis-grid{display:grid;grid-template-columns:repeat(3,1fr);gap:8px;margin-top:14px}.axis-grid article{border:1px solid #dbe2ea;border-radius:12px;padding:12px}.axis-grid small{display:block;color:#0f766e;font-weight:800;margin-top:4px}.independent{border:1px dashed #94a3b8;border-radius:12px;padding:12px}.candidate-value{font-weight:800;color:#0e7490}.boundary{border-left:4px solid #f59e0b;padding:9px 12px;background:#fffbeb}.priority{display:grid;grid-template-columns:32px 1fr;gap:10px;background:#fff;border:1px solid #dbe2ea;border-left:4px solid;border-radius:12px;padding:14px;margin:8px 0}.priority-number{width:26px;height:26px;border-radius:50%;background:#e2e8f0;display:grid;place-items:center;font-weight:800}.scroll{overflow-x:auto}table{width:100%;border-collapse:collapse;background:#fff}th,td{text-align:left;vertical-align:top;padding:9px;border-bottom:1px solid #e2e8f0}.mismatch{background:#fff1f2;color:#9f1239}dl{display:grid;grid-template-columns:minmax(150px,220px) 1fr;gap:6px 14px}dt{font-weight:700}dd{margin:0;font-family:monospace;overflow-wrap:anywhere}pre{background:#0f172a;color:#a7f3d0;padding:16px;border-radius:12px;overflow:auto;font:12px/1.6 monospace}ul{padding-left:20px}[data-theme="dark"]{color-scheme:dark}[data-theme="dark"] body{background:#0b0f19;color:#f1f5f9}[data-theme="dark"] h2{color:#94a3b8}[data-theme="dark"] .eyebrow{color:#34d399}[data-theme="dark"] .muted,[data-theme="dark"] small{color:#94a3b8}[data-theme="dark"] .outcome,[data-theme="dark"] .section-card{background:#1e293b;border-color:#334155}[data-theme="dark"] .callout{background:#451a03;border-color:#9a3412;color:#fed7aa}[data-theme="dark"] .axis-grid article{border-color:#475569}[data-theme="dark"] .axis-grid small{color:#5eead4}[data-theme="dark"] .candidate-value{color:#67e8f9}[data-theme="dark"] .boundary{background:#422006}[data-theme="dark"] .priority{background:#1e293b;border-color:#334155}[data-theme="dark"] .priority-number{background:#334155;color:#f8fafc}[data-theme="dark"] table{background:#1e293b}[data-theme="dark"] th,[data-theme="dark"] td{border-bottom:1px solid #334155}[data-theme="dark"] .mismatch{background:#4c0519;color:#fda4af}[data-theme="dark"] pre{background:#020617;color:#6ee7b7;border:1px solid #1e293b}@media(max-width:640px){.page{padding:24px 14px}dl{grid-template-columns:1fr}.axis-grid{grid-template-columns:1fr}dd{margin-bottom:7px}}@media print{body{background:#fff!important;color:#0f172a!important}.page{padding:0}.outcome,.section-card,.priority,table{background:#fff!important;break-inside:avoid}.callout{border-color:#999}}
</style></head><body><main class="page"><p class="eyebrow">FactorForge CDS Design Review</p><h1>Design Comparison &amp; Decision Report</h1><p class="muted">Result ${escapeHtml(reportValue(model.identity.result_id))} · created ${escapeHtml(reportValue(model.identity.result_created_at))} · report generated ${escapeHtml(model.identity.report_generated_at)}</p><div class="callout"><b>Sequence-data notice:</b> This HTML contains ${referenceSequence ? 'the reference and optimized DNA sequences' : 'the optimized DNA sequence'}. Handle and share it according to your sequence-data policy.</div><section><h2>Decision brief</h2><div class="outcome" style="border-color:${decisionColor}"><strong style="color:${decisionColor}">${escapeHtml(model.disposition.automated_decision.replaceAll('_', ' '))}</strong><p>${escapeHtml(model.interpretation.headline)}</p><p class="counts">${escapeHtml(`${reportValue(model.disposition.required_failure_count)} required fail · ${reportValue(model.disposition.preferred_warning_count)} warning · ${reportValue(model.disposition.unavailable_check_count)} unavailable`)}</p><small>${escapeHtml(reportValue(model.disposition.explanation))}</small></div></section>${structure3dReport.html}${standaloneContract}${scientificReviewHtml(model)}${standaloneComparison}<section><h2>Review priorities and next actions</h2>${priorities}</section><section><h2>Requested vs applied settings</h2><div class="scroll"><table><thead><tr><th>Setting</th><th>Requested</th><th>Applied / recorded</th></tr></thead><tbody>${settings}</tbody></table></div></section><section><h2>All computational checks</h2><div class="scroll"><table><thead><tr><th>Detection</th><th>Enforcement</th><th>Action</th><th>Status</th></tr></thead><tbody>${checks}</tbody></table></div></section>${candidates}<section><h2>Sequence and process</h2><p>Input type: <b>${escapeHtml(model.context.input_type)}</b> · output length: <b>${escapeHtml(reportValue(model.sequence_summary.output_length_nt, ' nt'))}</b> · nucleotide comparison: <b>${escapeHtml(model.sequence_summary.comparison_available ? reportValue(model.sequence_summary.nucleotide_changes) : 'Not recorded')}</b></p><p>Domestication: ${escapeHtml(model.process.domestication_attempted ? 'Attempted' : 'Not attempted')} · MFE: ${escapeHtml(model.metrics.mfe_status === 'computed' ? 'Computed' : `Not computed (${model.metrics.mfe_status_reason})`)}</p></section><section><h2>Downstream handoff checklist</h2><ul>${model.interpretation.next_steps.map(item => `<li>${escapeHtml(item)}</li>`).join('')}</ul></section><section><h2>Reproducibility and provenance</h2><dl>${provenanceRows}</dl></section>${referenceSequence ? `<section><h2>Reference sequence (DNA)</h2><pre>${escapeHtml(referenceSequence)}</pre></section>` : ''}<section><h2>Optimized sequence (DNA)</h2><pre>${escapeHtml(sequence)}</pre></section><section><h2>Interpretation and limitations</h2><p>${escapeHtml(model.interpretation.scope)}</p><ul>${model.interpretation.limitations.map(item => `<li>${escapeHtml(item)}</li>`).join('')}</ul></section></main>${structure3dReport.script}</body></html>`;
}

function downloadResultsReportHtml(model) {
    trackEvent('report_download', { format: 'html' });
    downloadTextArtifact(standaloneReportHtml(model), 'text/html;charset=utf-8', `factorforge_design_review_${reportFileStem(model)}.html`);
    showToast('Design review report downloaded', 'success');
}

function standaloneInteractiveReportHtml(model) {
    const accounting = model.accounting || { available: false, codons: [] };
    const safeData = JSON.stringify({ accounting }).replace(/</g, '\\u003c');
    const summary = accounting.available
        ? `<div class="accounting"><article><small>Total sense codons</small><b>${accounting.totalCodons}</b></article><article><small>Synonymous shifts</small><b>${accounting.synonymousChanges} · ${percentOf(accounting.synonymousChanges, accounting.totalCodons)}</b></article><article><small>Unchanged</small><b>${accounting.unchanged} · ${percentOf(accounting.unchanged, accounting.totalCodons)}</b></article><article><small>AA preservation</small><b>${accounting.aaIdentity.toFixed(1)}% · ${accounting.substitutions} substitutions</b></article></div>`
        : '<p class="callout">Codon-by-codon accounting is unavailable because the run did not include a comparable input CDS.</p>';
    const interactive = `<section><p class="eyebrow">Same Protein. Better-Designed DNA.</p><h2>Codon shift accounting and interactive map</h2><p>FactorForge redesigns synonymous codons and sequence-level features while verifying the encoded amino-acid sequence. Values below are calculated from the embedded input and output CDS; unavailable host frequencies remain N/A.</p>${summary}<div class="standalone-track"><div><span>Original</span><canvas id="reportOriginal" height="56"></canvas></div><div><span>Optimized</span><canvas id="reportOptimized" height="56"></canvas></div><div id="reportInspector" class="report-inspector">${accounting.available ? 'Move across either track to inspect a codon.' : 'No comparable CDS track.'}</div></div></section>`;
    const script = `<script>const reportData=${safeData};(()=>{const a=reportData.accounting,cs=a.codons||[],o=document.getElementById('reportOriginal'),p=document.getElementById('reportOptimized'),i=document.getElementById('reportInspector');const color=f=>f==null?'#334155':f*100<20?'#ef4444':f*100<35?'#f59e0b':f*100<50?'#64748b':f*100<=70?'#22c55e':'#10b981';const fmt=f=>f==null?'N/A':(f*100).toFixed(1)+'%';function draw(c,key,h){const r=c.getBoundingClientRect(),d=devicePixelRatio||1,w=Math.max(1,r.width);c.width=w*d;c.height=56*d;const x=c.getContext('2d');x.scale(d,d);x.fillStyle='#020617';x.fillRect(0,0,w,56);if(!cs.length)return;const sw=w/cs.length;cs.forEach((v,n)=>{x.fillStyle=color(v[key]);x.fillRect(n*sw,13,Math.max(1,sw+.25),30)});if(h!=null){x.strokeStyle='#fff';x.beginPath();x.moveTo((h+.5)*sw,2);x.lineTo((h+.5)*sw,54);x.stroke()}}function render(h){draw(o,'origFreq',h);draw(p,'optFreq',h)}function move(e){if(!cs.length)return;const r=p.getBoundingClientRect(),n=Math.max(0,Math.min(cs.length-1,Math.floor((e.clientX-r.left)/r.width*cs.length))),v=cs[n];render(n);i.textContent='Codon '+(n+1)+' · nt '+(n*3+1)+'–'+(n*3+3)+' · '+v.origCodon+' ('+fmt(v.origFreq)+') → '+v.optCodon+' ('+fmt(v.optFreq)+') · AA '+v.aa+' · local GC '+(v.localGc==null?'N/A':v.localGc.toFixed(1)+'%')+' · '+(v.isChanged?'change observed':'unchanged')}[o,p].forEach(c=>c.addEventListener('mousemove',move));addEventListener('resize',()=>render(null));render(null)})();<\/script>`;
    return standaloneReportHtml(model)
        .replace('</style>', '.accounting{display:grid;grid-template-columns:repeat(4,1fr);gap:8px}.accounting article{background:#fff;border:1px solid #dbe2ea;border-radius:12px;padding:12px}.accounting b,.accounting small{display:block}.standalone-track{background:#0f172a;border-radius:14px;padding:16px}.standalone-track>div{display:grid;grid-template-columns:70px 1fr;align-items:center;gap:8px;margin:8px 0;color:#cbd5e1;font-size:12px;font-weight:700}.standalone-track canvas{width:100%;height:56px;border-radius:8px;cursor:crosshair}.report-inspector{grid-template-columns:1fr!important;background:#1e293b;padding:12px;border-radius:8px}@media(max-width:640px){.accounting{grid-template-columns:1fr 1fr}}</style>')
        .replace('</main>', `${interactive}</main>`)
        .replace('</body>', `${script}</body>`);
}

function downloadInteractiveReport(model) {
    trackEvent('report_download', { format: 'interactive_html' });
    const constructId = model.identity.construct_id || model.identity.result_id || 'Construct';
    const safeId = String(constructId).replace(/[^A-Za-z0-9._-]+/g, '_');
    downloadTextArtifact(standaloneInteractiveReportHtml(model), 'text/html;charset=utf-8', `FactorForge_Optimization_Report_${safeId}.html`);
    showToast('Interactive report downloaded', 'success');
}

async function downloadEvidencePackageZip() {
    if (!state.results) return showToast('Run an optimization before downloading the evidence package.', 'info');
    if (typeof JSZip === 'undefined') {
        return showToast('ZIP compiler library loading, please try again.', 'error');
    }
    trackEvent('evidence_package_download', { target: state.sequenceName || 'target' });

    try {
        showToast('Compiling complete evidence package (.zip)...', 'info');
        const zip = new JSZip();
        const fastaFolder = zip.folder('sequences');
        const primary = getPrimaryResult(state.results);
        const optSeq = primary.optimized_sequence;
        const targetName = (state.sequenceName || 'Target_Sequence').replace(/[^a-zA-Z0-9_-]/g, '_');
        const timestamp = new Date().toISOString();

        // 1. FASTA files
        const policyLabel = state.codonPolicy === 'plantform_balanced_codon_v1' ? 'balanced_non_cpg' : (state.codonPolicy || 'standard');
        const candidates = (state.results.candidates || []).filter(candidate => candidate.dna_sequence);
        const multiFasta = candidates.length ? candidates.map(candidate => `>${targetName}_${String(candidate.id).replace(/[^A-Za-z0-9_-]/g,'_')} | Policy=${policyLabel}\n${candidate.dna_sequence}\n`).join('') : `>${targetName}_OPTIMIZED | GC=${primary.metrics.gc_percent ?? 'N/A'}% | Policy=${policyLabel}\n${optSeq}\n`;
        fastaFolder.file("all_variants.fasta", multiFasta);
        fastaFolder.file(`${targetName}_optimized.fasta`, `>${targetName}_OPTIMIZED | Policy=${policyLabel}\n${optSeq}\n`);
        const reference = state.reportData?.artifacts?.reference_sequence;
        if (reference) {
            fastaFolder.file(`${targetName}_original.fasta`, `>${targetName}_ORIGINAL\n${reference}\n`);
        }

        // 2. Structured JSON
        const expJson = {
            experiment_id: `EXP-${Date.now()}-${targetName}`,
            target_name: targetName,
            host_organism: state.host || 'nbenthamiana',
            timestamp: timestamp,
            results: state.results,
            report_data: state.reportData || null,
        };
        zip.file("experiment.json", JSON.stringify(expJson, null, 2));

        // 3. Markdown Design Dossier
        const mdDossier = evidenceDossierMarkdown(state.results, primary, targetName, timestamp);
        zip.file("DESIGN_DOSSIER.md", mdDossier);

        // 4. Standalone HTML Report
        if (state.reportData) {
            const reportHtml = standaloneInteractiveReportHtml(state.reportData);
            zip.file("dashboard.html", reportHtml);
        }

        // 5. Generate and download ZIP
        const content = await zip.generateAsync({ type: "blob" });
        const url = window.URL.createObjectURL(content);
        const a = document.createElement('a');
        a.href = url;
        a.download = `FactorForge_Evidence_Package_${targetName}_${Date.now()}.zip`;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        window.URL.revokeObjectURL(url);
        showToast('Evidence package (.zip) downloaded successfully!', 'success');
    } catch (err) {
        console.error('Evidence package compilation error:', err);
        showToast('Failed to compile evidence package: ' + err.message, 'error');
    }
}

function evidenceDossierMarkdown(result, primary, targetName, timestamp) {
    const typeIis = (result.qc_decision_matrix || []).find(check => check.criterion === 'type_iis');
    const observation = result.scientific_review || {};
    return [
        `# FactorForge Design Dossier: ${targetName}`,
        ``,
        `- **Date:** ${timestamp}`,
        `- **Host Organism:** ${result.host_profile || 'Not recorded'}`,
        `- **Product Version:** ${result.product_version || 'Not recorded'}`,
        `- **Engine:** ${result.cds_design?.engine || 'Not recorded'} · ${result.cds_design?.engine_version || 'Not recorded'}`,
        ``,
        `## 1. Sequence Metrics Summary`,
        ``,
        `| Metric | Value |`,
        `|---|---|`,
        `| Sequence Length | ${primary.optimized_sequence.length} nt |`,
        `| GC Content | ${primary.metrics.gc_percent ?? 'N/A'}% |`,
        `| CAI Score | ${primary.metrics.cai ?? 'N/A'} |`,
        `| Recorded amino-acid identity | ${result.constraint_report?.aa_identity ?? 'Not recorded'} |`,
        ``,
        `## 2. Computational Review Summary`,
        `- Recorded Type IIS check: ${typeIis?.result || 'Not recorded'}`,
        `- Policy decision: ${result.policy_decision || 'Not recorded'}`,
        `- RNA scan: ${observation.rna_risks?.summary_status || 'Not recorded'}`,
        `- Initiation-window MFE: ${observation.initiation_mfe?.status || 'Not recorded'}`,
        `- Computational observations do not establish laboratory acceptance, expression or folding benefits.`,
    ].join('\n');
}

function downloadEvidenceRecordJson(model) {
    trackEvent('report_download', { format: 'evidence_json' });
    const content = `${JSON.stringify(evidenceRecordFromModel(model), null, 2)}\n`;
    downloadTextArtifact(content, 'application/json;charset=utf-8', `factorforge_design_evidence_${reportFileStem(model)}.json`);
    showToast('Sequence-free evidence record downloaded', 'success');
}

function renderCustomRestrictionResults(res) {
    if (!elements.customRestrictionResults || !elements.customRestrictionResultsBody) return;

    const custom = res.custom_restriction_sites;
    if (!custom) {
        elements.customRestrictionResults.classList.remove('hidden', 'domestication-attempted');
        elements.customRestrictionResults.classList.add('domestication-not-attempted');
        elements.customRestrictionResultsBody.innerHTML = '<p class="font-bold text-slate-600 dark:text-slate-300">Domestication not attempted</p>';
        return;
    }

    elements.customRestrictionResults.classList.remove('domestication-not-attempted');
    elements.customRestrictionResults.classList.add('domestication-attempted');

    const removed = Array.isArray(custom.removed) ? custom.removed : [];
    const unresolved = Array.isArray(custom.unresolved) ? custom.unresolved : [];
    const before = res.metrics && res.metrics.before;
    const after = res.metrics && res.metrics.after;

    const removedHtml = removed.length > 0
        ? removed.map(site => `
            <li class="flex items-start justify-between gap-3 py-1">
                <span><span class="text-emerald-600 font-black">✓</span> ${escapeHtml(site.name)} at ${site.position}</span>
                <span class="font-mono text-[11px] text-emerald-700 dark:text-emerald-300">${escapeHtml(site.substitution || '')}</span>
            </li>
        `).join('')
        : '<li class="py-1 text-slate-500 dark:text-slate-400">No custom sites removed</li>';

    const unresolvedHtml = unresolved.length > 0
        ? unresolved.map(site => `
            <li class="flex items-start justify-between gap-3 py-1">
                <span><span class="text-amber-600 font-black">!</span> ${escapeHtml(site.name)} at ${site.position}</span>
                <span class="text-[11px] text-amber-700 dark:text-amber-300">${escapeHtml(site.reason || 'unresolved')}</span>
            </li>
        `).join('')
        : '<li class="py-1 text-slate-500 dark:text-slate-400">No unresolved custom sites</li>';

    const metricsHtml = before && after ? `
        <div class="grid grid-cols-2 gap-3 pt-3 border-t border-slate-100 dark:border-slate-800 text-xs">
            <div class="rounded-xl bg-slate-50 dark:bg-slate-800 p-3">
                <p class="text-[10px] font-extrabold uppercase tracking-widest text-slate-500">Before</p>
                <p class="mt-1 font-bold text-slate-800 dark:text-slate-100">CAI ${Number(before.cai || 0).toFixed(3)} · GC ${Number(before.gc || 0).toFixed(1)}%</p>
            </div>
            <div class="rounded-xl bg-emerald-50 dark:bg-emerald-900/20 p-3">
                <p class="text-[10px] font-extrabold uppercase tracking-widest text-emerald-600">After</p>
                <p class="mt-1 font-bold text-slate-800 dark:text-slate-100">CAI ${Number(after.cai || 0).toFixed(3)} · GC ${Number(after.gc || 0).toFixed(1)}%</p>
            </div>
        </div>
    ` : '';

    elements.customRestrictionResultsBody.innerHTML = `
        <div class="grid grid-cols-1 gap-4 text-xs">
            <p class="font-bold text-emerald-700 dark:text-emerald-300">Domestication attempted</p>
            <div>
                <p class="text-[10px] font-extrabold uppercase tracking-widest text-emerald-600 mb-2">Removed sites</p>
                <ul class="divide-y divide-slate-100 dark:divide-slate-800">${removedHtml}</ul>
            </div>
            <div>
                <p class="text-[10px] font-extrabold uppercase tracking-widest text-amber-600 mb-2">Unresolved sites</p>
                <ul class="divide-y divide-slate-100 dark:divide-slate-800">${unresolvedHtml}</ul>
            </div>
            ${metricsHtml}
        </div>
    `;
    elements.customRestrictionResults.classList.remove('hidden');
}

function escapeHtml(value) {
    return String(value)
        .replace(/&/g, '&amp;')
        .replace(/</g, '&lt;')
        .replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;')
        .replace(/'/g, '&#039;');
}


function updateValidationIcon(id, pass) {
    const el = document.getElementById(id);
    if (el) {
        if (pass === null) {
            el.textContent = '⚠️';
            el.className = 'text-amber-500';
        } else {
            el.textContent = pass ? '✅' : '❌';
            el.className = pass ? 'text-emerald-400' : 'text-rose-500';
        }
    }
}

const DOMAIN_GROUP_LABELS = {
    configured_constraint: 'Configured constraint',
    advisory_scan: 'Advisory sequence scans',
    assembly_review: 'Assembly review',
};

function statusIcon(status) {
    if (status === 'PASS') return { icon: '✅', className: 'text-emerald-400' };
    if (status === 'NOT_RUN' || status === 'NOT_APPLICABLE') return { icon: '⏭️', className: 'text-slate-500' };
    if (status === 'WARNING') return { icon: '⚠️', className: 'text-amber-500' };
    return { icon: '❌', className: 'text-rose-500' };
}

function resultChecksFromResponse(res) {
    if (res.validation && res.validation.checks) return res.validation.checks;
    if (res.validation_report && res.validation_report.checks) return res.validation_report.checks;
    return {};
}

function renderValidationChecks(registryChecks, resultChecks) {
    const container = document.getElementById('validationChecksContainer');
    if (!container) return;

    if (!Array.isArray(registryChecks) || registryChecks.length === 0) {
        container.innerHTML = '<p class="text-xs text-slate-500">Validation registry unavailable.</p>';
        return;
    }

    const sorted = [...registryChecks].sort((a, b) => a.order - b.order);
    const groups = new Map();
    for (const check of sorted) {
        const groupLabel = DOMAIN_GROUP_LABELS[check.primary_domain] || check.primary_domain;
        if (!groups.has(groupLabel)) groups.set(groupLabel, []);
        groups.get(groupLabel).push(check);
    }

    let html = '';
    for (const [groupLabel, checks] of groups) {
        html += `<div>
            <p class="text-[9px] font-bold text-slate-600 uppercase tracking-widest mb-2">${escapeHtml(groupLabel)}</p>
            <div class="space-y-2">`;
        for (const check of checks) {
            const result = resultChecks[check.check_id];
            const status = result ? result.status : 'NOT_RUN';
            const { icon, className } = statusIcon(status);
            const findingText = result && result.finding_count != null ? ` (${result.finding_count})` : '';
            html += `<div class="flex items-center space-x-3 text-xs text-slate-300">
                <span class="${className}">${icon}</span>
                <span class="font-medium">${escapeHtml(check.display_name)}${findingText}</span>
            </div>`;
        }
        html += `</div></div>`;
    }

    container.innerHTML = html;
}

function formatSequence(seq) {
    if (!seq) return '';
    return seq.match(/.{1,60}/g).join('\n');
}

function renderGCGraph(seq, hostId = state.host) {
    const windowSize = 50;
    const data = [];
    const labels = [];

    for (let i = 0; i < seq.length - windowSize; i += 10) {
        const window = seq.substring(i, i + windowSize);
        const gcCount = (window.match(/[GC]/gi) || []).length;
        data.push((gcCount / windowSize) * 100);
        labels.push(i);
    }

    if (chartInstance) chartInstance.destroy();

    const isDark = document.documentElement.classList.contains('dark');
    const textColor = isDark ? '#94a3b8' : '#475569';
    const gridColor = isDark ? 'rgba(255,255,255,0.08)' : 'rgba(0,0,0,0.05)';
    const n = labels.length;
    // Host-aware reference band (v3.3.0) — previously hardcoded to
    // 55/65 for every host, which mislabeled the N. benthamiana v2 default.
    const { gc_min: bandMin, gc_max: bandMax } = getGcRange(hostId);
    if (elements.gcZoneLabel) {
        elements.gcZoneLabel.textContent = `Reference Band ${bandMin}–${bandMax}%`;
    }

    const ctx = elements.gcChart.getContext('2d');
    chartInstance = new Chart(ctx, {
        type: 'line',
        data: {
            labels: labels,
            datasets: [
                {
                    label: `Target Max (${bandMax}%)`,
                    data: Array(n).fill(bandMax),
                    borderColor: 'rgba(16, 185, 129, 0.35)',
                    borderWidth: 1,
                    borderDash: [4, 4],
                    pointRadius: 0,
                    fill: false,
                    tension: 0
                },
                {
                    label: `Target Min (${bandMin}%)`,
                    data: Array(n).fill(bandMin),
                    borderColor: 'rgba(16, 185, 129, 0.35)',
                    borderWidth: 1,
                    borderDash: [4, 4],
                    pointRadius: 0,
                    fill: '-1',
                    backgroundColor: 'rgba(16, 185, 129, 0.06)',
                    tension: 0
                },
                {
                    label: 'GC Content %',
                    data: data,
                    borderColor: '#10B981',
                    backgroundColor: 'rgba(16, 185, 129, 0.15)',
                    borderWidth: 2,
                    pointRadius: 0,
                    fill: false,
                    tension: 0.4
                }
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                legend: { display: false },
                tooltip: {
                    mode: 'index',
                    intersect: false,
                    filter: item => item.datasetIndex === 2
                }
            },
            scales: {
                y: {
                    min: 20,
                    max: 80,
                    grid: { color: gridColor },
                    ticks: { color: textColor, font: { size: 10 } }
                },
                x: {
                    display: false
                }
            }
        }
    });
}

function historyResultSnapshot(result, primary) {
    const clone = value => value == null ? value : JSON.parse(JSON.stringify(value));
    const candidates = Array.isArray(result.candidates) ? result.candidates.map(candidate => ({
        id: candidate.id || null,
        label: candidate.label || candidate.id || null,
        cai: candidate.cai ?? null,
        gc_percent: candidate.gc_percent ?? null,
        dna_sequence: candidate.dna_sequence || null,
        automated_decision: candidate.automated_decision || null,
        required_failure_count: candidate.required_failure_count ?? null,
        preferred_warning_count: candidate.preferred_warning_count ?? null,
    })) : [];
    const lightViewer = result.structure_viewer ? {
        source_type: result.structure_viewer.source_type || null,
        source_label: result.structure_viewer.source_label || null,
        positions: result.structure_viewer.positions || [],
        ca_trace: result.structure_viewer.ca_trace || [],
    } : null;
    return {
        optimized_sequence: primary.optimized_sequence,
        original_length: result.original_length ?? result.input_summary?.length ?? null,
        optimized_length: result.optimized_length ?? primary.optimized_sequence?.length ?? null,
        metrics: clone(result.metrics || primary.metrics || {}),
        validation: clone(result.validation || primary.validation || {}),
        profile: result.profile || state.objective,
        host_profile: getResultHostProfile(result),
        input_type: result.input_type || result.validation?.input_type || result.provenance?.normalized_input_type || null,
        input_summary: clone(result.input_summary || null),
        acceptance_criteria_snapshot: clone(result.acceptance_criteria_snapshot || null),
        automated_decision: result.automated_decision || null,
        decision_summary: clone(result.decision_summary || null),
        qc_decision_matrix: clone(result.qc_decision_matrix || null),
        acceptance_evaluation: clone(result.acceptance_evaluation || null),
        reviewer_disposition: clone(result.reviewer_disposition || null),
        result_identifier: result.result_identifier || null,
        construct_id: result.construct_id || null,
        created_at: result.created_at || null,
        product_version: result.product_version || null,
        reference_policy_version: result.reference_policy_version || null,
        codon_reference_id: result.codon_reference_id || null,
        gc_reference_band: result.gc_reference_band || null,
        metadata: clone(result.metadata || null),
        provenance: clone(result.provenance || null),
        cds_design: clone(result.cds_design || null),
        constraint_report: clone(result.constraint_report || null),
        custom_restriction_sites: clone(result.custom_restriction_sites || null),
        scientific_review: clone(result.scientific_review || null),
        candidate_scientific_reviews: Array.isArray(result.candidates) ? result.candidates.map(c=>({id:c.id,label:c.label,scientific_review:clone(c.scientific_review || null)})) : clone(result.candidate_scientific_reviews || []),
        design_mode: result.design_mode || null,
        structure_comparison: clone(result.structure_comparison || null),
        structure_viewer: lightViewer,
        candidates: clone(result.candidates || null),
        report_candidates: candidates,
    };
}

function addToHistory(input, result) {
    const primary = getPrimaryResult(result);
    const item = {
        id: Date.now(),
        schemaVersion: HISTORY_SCHEMA_VERSION,
        timestamp: new Date().toLocaleString(),
        inputLen: input.length,
        inputType: result.input_type || result.validation?.input_type || result.provenance?.normalized_input_type || null,
        profile: result.design_mode === 'paired_3dp' ? 'structure_3dp' : state.objective,
        host: getResultHostProfile(result),
        cai: primary.metrics.cai,
        gc: calculateGC(primary.optimized_sequence),
        sequence: primary.optimized_sequence,
        resultSnapshot: historyResultSnapshot(result, primary),
    };

    state.history.unshift(item);
    if (state.history.length > 10) state.history.pop();
    try {
        localStorage.setItem('factorforge_history', JSON.stringify({ schemaVersion: HISTORY_SCHEMA_VERSION, items: state.history }));
    } catch (err) {
        console.warn('Unable to persist full history to localStorage:', err);
    }
    renderHistory();
}

function renderHistory() {
    if (!elements.historyList) return;

    if (state.history.length === 0) {
        elements.historyList.innerHTML = '<p class="text-[10px] text-slate-400 text-center py-4">No recent history</p>';
        return;
    }

    elements.historyList.innerHTML = state.history.map(item => `
        <div class="p-2 border border-slate-100 rounded-lg hover:bg-slate-50 cursor-pointer transition-all mb-2 flex justify-between items-center group" onclick="loadHistoryItem(${item.id})">
            <div>
                <p class="text-[10px] font-bold text-slate-700">${item.inputLen}${item.inputType === 'protein' ? 'aa' : 'bp'} → ${item.profile}</p>
                <p class="text-[9px] text-slate-400">${formatHostProfile(item.host || 'nbenthamiana')} · ${item.timestamp}</p>
            </div>
            <div class="text-right">
                <p class="text-[10px] font-bold text-emerald-600">CAI: ${item.cai.toFixed(3)}</p>
                <p class="text-[9px] text-slate-400">GC: ${item.gc}%</p>
            </div>
        </div>
    `).join('');
}

window.loadHistoryItem = (id) => {
    const item = state.history.find(h => h.id === id);
    if (!item) return;

    if (item.resultSnapshot) {
        elements.sequenceInput.value = '';
        state.sequence = '';
        elements.previewContainer.classList.add('hidden');
        elements.inputTypeBadge.classList.add('hidden');
        updateInputStats('');
        state.results = JSON.parse(JSON.stringify(item.resultSnapshot));
    } else {
        const legacyInput = item.inputSequence || '';
        elements.sequenceInput.value = legacyInput;
        handleSequenceChange({ target: { value: legacyInput } });
        state.results = {
            optimized_sequence: item.sequence,
            metrics: { cai: item.cai, gc_percent: item.gc, polya_signals: 0, length: item.sequence.length },
            profile: item.profile,
            host_profile: item.host || 'nbenthamiana',
            acceptance_criteria_snapshot: item.acceptanceCriteria || null,
            automated_decision: item.automatedDecision || null,
            reviewer_disposition: item.reviewerDisposition || null,
        };
    }
    renderResults();
    showToast(item.resultSnapshot ? 'History report loaded (input sequence not stored)' : 'Legacy history item loaded', 'success');
};

function clearHistory() {
    state.history = [];
    localStorage.removeItem('factorforge_history');
    renderHistory();
    showToast('History cleared', 'info');
}

function clearAll() {
    invalidateStructureSuggestion();
    elements.sequenceInput.value = '';
    elements.fileUpload.value = '';
    elements.customRestrictionSites.value = '';
    elements.optimizationSeed.value = '';
    Array.from(elements.typeIisEnzymes).forEach(input => { input.checked = false; });
    state.sequence = '';
    state.customRestrictionSites = [];
    state.selectedTypeIisEnzymes = [];
    state.results = null;
    elements.previewContainer.classList.add('hidden');
    elements.inputTypeBadge.classList.add('hidden');
    updateInputStats('');
    elements.resultsContainer.classList.add('hidden');
    elements.constructIdDisplay.textContent = '';
    elements.constructIdRow.classList.add('hidden');
    if (elements.candidateComparisonContainer) elements.candidateComparisonContainer.classList.add('hidden');
    if (elements.customRestrictionResults) elements.customRestrictionResults.classList.add('hidden');
    if (elements.mfeWarningBanner) elements.mfeWarningBanner.classList.add('hidden');
    elements.emptyState.classList.remove('hidden');
    elements.validationStatus.classList.add('hidden');
    elements.canvasOriginalTrack?.classList.remove('has-track');
    elements.canvasOptimizedTrack?.classList.remove('has-track');
    updateDesignBriefSummary();

    showToast('Input cleared', 'info');
}

function toggleDetailsPanel() {
    const isHidden = elements.detailsContent.classList.contains('hidden');
    if (isHidden) {
        elements.detailsContent.classList.remove('hidden');
        elements.toggleArrow.style.transform = 'rotate(90deg)';
    } else {
        elements.detailsContent.classList.add('hidden');
        elements.toggleArrow.style.transform = 'rotate(0deg)';
    }
}

function toggleChangelog() {
    const isHidden = elements.changelogModal.classList.toggle('hidden');
    if (!isHidden) {
        const scrollArea = elements.changelogModal.querySelector('.overflow-y-auto');
        if (scrollArea) scrollArea.scrollTop = 0;
    }
}

// Utilities
function debounce(func, wait) {
    let timeout;
    return function executedFunction(...args) {
        const later = () => {
            clearTimeout(timeout);
            func(...args);
        };
        clearTimeout(timeout);
        timeout = setTimeout(later, wait);
    };
}

function showToast(message, type = 'info') {
    const toast = document.createElement('div');
    toast.className = `toast-enter p-4 rounded-2xl shadow-2xl border flex items-center space-x-3 transition-all ${type === 'success' ? 'bg-emerald-50 border-emerald-200 text-emerald-800' :
        type === 'error' ? 'bg-rose-50 border-rose-200 text-rose-800' :
            'bg-blue-50 border-blue-200 text-blue-800'
        }`;

    const icon = type === 'success' ? '✅' : type === 'error' ? '🚫' : 'ℹ️';

    toast.innerHTML = `
        <span class="text-xl">${icon}</span>
        <span class="text-xs font-bold uppercase tracking-tight">${message}</span>
    `;

    elements.toastContainer.appendChild(toast);

    setTimeout(() => {
        toast.classList.replace('toast-enter', 'toast-exit');
        setTimeout(() => toast.remove(), 300);
    }, 4000);
}

async function copyToClipboard() {
    if (!state.results) return;
    try {
        await navigator.clipboard.writeText(getPrimaryResult(state.results).optimized_sequence);
        const originalText = elements.copyBtn.innerHTML;
        elements.copyBtn.innerHTML = '<span>🎉</span> <span>Copied!</span>';
        elements.copyBtn.classList.add('bg-emerald-100');
        setTimeout(() => {
            elements.copyBtn.innerHTML = originalText;
            elements.copyBtn.classList.remove('bg-emerald-100');
        }, 2000);
        showToast('Sequence copied to clipboard', 'success');
    } catch (err) {
        showToast('Failed to copy', 'error');
    }
}

async function copyConstructId() {
    const constructId = state.results?.construct_id;
    if (!constructId) return;
    try {
        await navigator.clipboard.writeText(constructId);
        const originalText = elements.copyConstructId.textContent;
        elements.copyConstructId.textContent = 'Copied';
        setTimeout(() => {
            elements.copyConstructId.textContent = originalText;
        }, 2000);
        showToast('Construct ID copied', 'success');
    } catch (err) {
        showToast('Failed to copy construct ID', 'error');
    }
}

function downloadFile(format) {
    if (!state.results) return;
    trackEvent('file_download', { format });

    let content = '';
    let fileName = '';
    const primary = getPrimaryResult(state.results);
    const seq = primary.optimized_sequence;
    const exportObjective = state.results.design_mode === 'paired_3dp' ? 'structure_3dp' : state.objective;

    if (format === 'fasta') {
        content = `>FactorForge_Optimized | Objective: ${exportObjective} | CAI: ${primary.metrics.cai}\n${seq}`;
        fileName = `optimized_sequence_${Date.now()}.fasta`;
    } else {
        // Basic GenBank template
        content = `LOCUS       Exported                ${seq.length} bp    DNA     linear   \n`;
        const hostLabel = HOST_LABELS[state.host] || 'N. benthamiana';
        content += `DEFINITION  FactorForge Optimized Sequence for ${hostLabel}\n`;
        content += `FEATURES             Location/Qualifiers\n`;
        content += `     CDS             1..${seq.length}\n`;
        content += `                     /label="Optimized_CDS"\n`;
        content += `                     /note="Objective: ${exportObjective}"\n`;
        content += `ORIGIN      \n`;

        const lines = seq.toLowerCase().match(/.{1,60}/g);
        lines.forEach((line, i) => {
            const start = (i * 60) + 1;
            const groups = line.match(/.{1,10}/g).join(' ');
            content += `${start.toString().padStart(9, ' ')} ${groups}\n`;
        });
        content += `//`;
        fileName = `optimized_sequence_${Date.now()}.gb`;
    }

    const blob = new Blob([content], { type: 'text/plain' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = fileName;
    a.click();
    window.URL.revokeObjectURL(url);
    showToast(`File ${fileName} ready`, 'success');
}

function submitValidation() {
    trackEvent('validation_submit', { objective: state.objective });
    const ISSUE_URL = 'https://github.com/eijex/factorforge-cds/issues/new';
    const params = new URLSearchParams({ template: 'wet_lab_result.yml' });

    if (state.results) {
        const version = state.results.engine_versions?.product || '3.5.7';
        const profile = state.results?.profile || state.objective || '';
        params.set('title', `[wet-lab-summary] ${version} ${profile}`.trim());
    }

    window.open(`${ISSUE_URL}?${params.toString()}`, '_blank', 'noopener,noreferrer');
}

async function copyJson() {
    const text = document.getElementById('jsonDetails').textContent;
    if (!text) return;
    await navigator.clipboard.writeText(text);
    const btn = elements.copyJsonBtn;
    btn.textContent = 'Copied!';
    btn.classList.add('text-emerald-400');
    setTimeout(() => { btn.textContent = 'Copy'; btn.classList.remove('text-emerald-400'); }, 2000);
}

// Static Data
function getMockResult() {
    const mockSeq = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCTTCAGCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCACGGCATGGACGAGCTGTACAAG";
    return {
        optimized_sequence: mockSeq,
        original_length: state.sequence.length,
        optimized_length: mockSeq.length,
        metrics: {
            cai: 0.884,
            gc_percent: 42.6,
            polya_signals: 0,
            length: mockSeq.length
        },
        profile: state.objective,
        host_profile: state.host,
        validation: {
            polya: 'PASS',
            moclo: 'UNCHECKED',
            gc: 'PASS'
        }
    };
}
function calculateGC(seq) {
    if (!seq) return 0;
    const gCount = (seq.match(/G/g) || []).length;
    const cCount = (seq.match(/C/g) || []).length;
    return parseFloat(((gCount + cCount) / seq.length * 100).toFixed(1));
}



// -----------------------------------------------------------------------------
// Upload-first policy workflow; existing advanced settings remain optional.
document.addEventListener('DOMContentLoaded', () => {
    document.getElementById('buildSopButton')?.addEventListener('click', () => {
        const settings = document.getElementById('manualSopOverrides');
        if (!settings) return;
        settings.open = true;
        settings.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
    });
});


// ========================================================
// Guided 2-State Workspace Controller & Tab Logic
// ========================================================

function updateWorkflowStepper(viewState, isRunning = false) {
    const b1 = elements.stepInputBadge;
    const b2 = elements.stepPolicyBadge;
    const b3 = elements.stepResultsBadge;
    if (!b1 || !b2 || !b3) return;

    const c1 = document.getElementById('stepInputCircle') || b1.querySelector('span:first-child');
    const c2 = document.getElementById('stepPolicyCircle') || b2.querySelector('span:first-child');
    const c3 = document.getElementById('stepResultsCircle') || b3.querySelector('span:first-child');

    const t1 = document.getElementById('stepInputText');
    const t2 = document.getElementById('stepPolicyText');
    const t3 = document.getElementById('stepResultsText');

    [b1, b2, b3].forEach(b => {
        b.className = 'flex items-center space-x-1.5 px-2.5 py-1 rounded-lg transition-all focus:outline-none cursor-pointer';
    });

    if (isRunning) {
        b1.className += ' bg-emerald-900/50 text-emerald-300 border border-emerald-500/40 opacity-80';
        if (c1) {
            c1.className = 'w-4 h-4 rounded-full bg-emerald-400 text-slate-900 text-[10px] font-black flex items-center justify-center';
            c1.textContent = '✓';
        }
        b2.className += ' bg-emerald-900/50 text-emerald-300 border border-emerald-500/40 opacity-80';
        if (c2) {
            c2.className = 'w-4 h-4 rounded-full bg-emerald-400 text-slate-900 text-[10px] font-black flex items-center justify-center';
            c2.textContent = '✓';
        }
        b3.className += ' bg-rose-600 text-white font-black shadow-lg ring-2 ring-rose-400 animate-pulse';
        if (c3) {
            c3.className = 'w-4 h-4 rounded-full bg-white text-rose-700 text-[10px] font-black flex items-center justify-center shadow-sm';
            c3.textContent = '⚡';
        }
        if (t3) t3.textContent = 'Optimizing...';
    } else if (viewState === 'review') {
        b1.className += ' bg-emerald-900/60 text-emerald-200 border border-emerald-500/30 hover:bg-emerald-800/60';
        if (c1) {
            c1.className = 'w-4 h-4 rounded-full bg-emerald-400 text-slate-900 text-[10px] font-black flex items-center justify-center shadow-sm';
            c1.textContent = '✓';
        }
        b2.className += ' bg-emerald-900/60 text-emerald-200 border border-emerald-500/30 hover:bg-emerald-800/60';
        if (c2) {
            c2.className = 'w-4 h-4 rounded-full bg-emerald-400 text-slate-900 text-[10px] font-black flex items-center justify-center shadow-sm';
            c2.textContent = '✓';
        }
        b3.className += ' bg-rose-600 text-white font-black shadow-lg ring-2 ring-rose-400 scale-105';
        if (c3) {
            c3.className = 'w-4 h-4 rounded-full bg-white text-rose-700 text-[10px] font-black flex items-center justify-center shadow-sm';
            c3.textContent = '3';
        }
        if (t3) t3.textContent = 'Results';
    } else {
        b1.className += ' bg-rose-600 text-white font-black shadow-md ring-2 ring-rose-400';
        if (c1) {
            c1.className = 'w-4 h-4 rounded-full bg-white text-rose-700 text-[10px] font-black flex items-center justify-center shadow-sm';
            c1.textContent = '1';
        }
        b2.className += ' bg-white/15 text-white/90 hover:bg-white/25 hover:text-white font-semibold';
        if (c2) {
            c2.className = 'w-4 h-4 rounded-full bg-emerald-300 text-slate-900 text-[10px] font-black flex items-center justify-center shadow-sm';
            c2.textContent = '2';
        }
        b3.className += ' text-white/40 hover:text-white/70 hover:bg-white/10';
        if (c3) {
            c3.className = 'w-4 h-4 rounded-full bg-white/20 text-white/60 text-[10px] font-bold flex items-center justify-center';
            c3.textContent = '3';
        }
        if (t3) t3.textContent = 'Results';
    }
}

function setWorkspaceState(viewState) {
    if (!elements.designWorkspace) return;
    elements.designWorkspace.setAttribute('data-view-state', viewState);

    updateWorkflowStepper(viewState, state.isOptimizing);

    if (viewState === 'review' && elements.reviewConfigSummary) {
        const hostName = elements.hostSelect ? elements.hostSelect.options[elements.hostSelect.selectedIndex].text.split('(')[0].trim() : 'N. benthamiana';
        const presetName = elements.sopPresetSelect ? elements.sopPresetSelect.options[elements.sopPresetSelect.selectedIndex].text.split('(')[0].trim() : 'FactorForge Default';
        elements.reviewConfigSummary.textContent = state.results?.design_mode === 'paired_3dp'
            ? `Host: ${hostName} · 3DP experimental + 1DP comparator · direct testing, approval/evidence deferred`
            : `Host: ${hostName} · Profile: ${presetName} · Standard Safety Baseline`;
    }
}

function initWorkspaceTabs() {
    const tabs = [
        { btn: elements.tabBtnSummary, id: 'summary' },
        { btn: elements.tabBtnComparison, id: 'comparison' },
        { btn: elements.tabBtnDifferences, id: 'differences' },
        { btn: elements.tabBtnTrace, id: 'trace' }
    ];

    function activateTab(tabId) {
        tabs.forEach(t => {
            if (t.btn) {
                if (t.id === tabId) {
                    t.btn.classList.add('active');
                    t.btn.setAttribute('aria-selected', 'true');
                } else {
                    t.btn.classList.remove('active');
                    t.btn.setAttribute('aria-selected', 'false');
                }
            }
        });

        // Toggle panes
        const summarySections = [
            elements.sameProteinBanner,
            elements.codonAccounting,
            elements.mfeWarningBanner,
            elements.designContractSummary,
            elements.automatedDecisionCard,
            elements.validationStatus,
            elements.qcDecisionMatrix,
            document.querySelector('#resultsContainer table.w-full'),
            elements.customRestrictionResults,
            document.querySelector('#resultsContainer .grid.grid-cols-2'),
            elements.resultsReport
        ];

        if (tabId === 'summary') {
            summarySections.forEach(s => s && s.classList.remove('hidden'));
            if (elements.tabPaneComparison) elements.tabPaneComparison.classList.add('hidden');
            if (elements.comparisonDashboard) elements.comparisonDashboard.classList.add('hidden');
            if (elements.tabPaneTrace) elements.tabPaneTrace.classList.add('hidden');
        } else if (tabId === 'comparison') {
            summarySections.forEach(s => s && s.classList.add('hidden'));
            if (elements.tabPaneComparison) elements.tabPaneComparison.classList.remove('hidden');
            if (elements.comparisonDashboard) elements.comparisonDashboard.classList.add('hidden');
            if (elements.tabPaneTrace) elements.tabPaneTrace.classList.add('hidden');
        } else if (tabId === 'differences') {
            summarySections.forEach(s => s && s.classList.add('hidden'));
            if (elements.tabPaneComparison) elements.tabPaneComparison.classList.add('hidden');
            if (elements.comparisonDashboard) elements.comparisonDashboard.classList.remove('hidden');
            if (elements.tabPaneTrace) elements.tabPaneTrace.classList.add('hidden');
        } else if (tabId === 'trace') {
            summarySections.forEach(s => s && s.classList.add('hidden'));
            if (elements.tabPaneComparison) elements.tabPaneComparison.classList.add('hidden');
            if (elements.comparisonDashboard) elements.comparisonDashboard.classList.add('hidden');
            if (elements.tabPaneTrace) elements.tabPaneTrace.classList.remove('hidden');
        }
    }

    tabs.forEach(t => {
        if (t.btn) {
            t.btn.addEventListener('click', () => activateTab(t.id));
        }
    });

    if (elements.btnEditSettings) {
        elements.btnEditSettings.addEventListener('click', () => {
            setWorkspaceState('configure');
            if (elements.sequenceInput) {
                elements.sequenceInput.scrollIntoView({ behavior: 'smooth', block: 'start' });
            }
        });
    }

    if (elements.stepInputBadge) {
        elements.stepInputBadge.addEventListener('click', () => {
            setWorkspaceState('configure');
            if (elements.sequenceInput) {
                elements.sequenceInput.focus();
                elements.sequenceInput.scrollIntoView({ behavior: 'smooth', block: 'center' });
            }
        });
    }
    if (elements.stepPolicyBadge) {
        elements.stepPolicyBadge.addEventListener('click', () => {
            setWorkspaceState('configure');
            if (elements.designBriefPanel) {
                elements.designBriefPanel.scrollIntoView({ behavior: 'smooth', block: 'start' });
            }
        });
    }
    if (elements.stepResultsBadge) {
        elements.stepResultsBadge.addEventListener('click', () => {
            if (state.lastResult) {
                setWorkspaceState('review');
                if (elements.resultsPanel) {
                    elements.resultsPanel.scrollIntoView({ behavior: 'smooth', block: 'start' });
                }
            }
        });
    }
}

document.addEventListener('DOMContentLoaded', () => {
    initWorkspaceTabs();
    setWorkspaceState('configure');
});
