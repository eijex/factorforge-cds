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
    ntabacum: 'Tobacco BY-2'
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
    ntabacum: { gc_min: 55.0, gc_max: 65.0 }
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
let apiCapabilities = {};

function getGcRange(hostId) {
    return hostGcRanges[hostId] || OFFLINE_GC_RANGES[hostId] || OFFLINE_GC_RANGES.nbenthamiana;
}

let validationRegistry = [];
const HISTORY_SCHEMA_VERSION = 3;

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
    objective: 'feasibility_best',
    host: 'nbenthamiana',
    saveDb: false,
    alignmentPage: 0,
    useTemplate: false,
    kozak: false,
    dinuc: false,
    customRestrictionSites: [],
    selectedTypeIisEnzymes: [],
    reviewerDisposition: null,
    results: null,
    isOptimizing: false,
    history: loadVersionedHistory()
};

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
    downloadFasta: document.getElementById('downloadFasta'),
    downloadGenbank: document.getElementById('downloadGenbank'),
    copyBtn: document.getElementById('copyBtn'),
    constructIdRow: document.getElementById('constructIdRow'),
    constructIdDisplay: document.getElementById('constructIdDisplay'),
    copyConstructId: document.getElementById('copyConstructId'),
    submitValidationBtn: document.getElementById('submitValidationBtn'),
    alphafoldLink: document.getElementById('alphafoldLink'),
    esmatlasFoldLink: document.getElementById('esmatlasFoldLink'),
    copyJsonBtn: document.getElementById('copyJsonBtn'),
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
    await loadApiMetadata();
    initEventListeners();
    updateDesignBriefSummary();
    renderHistory();
    console.log('FactorForge v3.4.5 Engaged');
});

// Loads server-owned GC ranges and validation labels. Supported hosts remain
// available to API/CLI clients, while the public web workflow stays focused on
// the production N. benthamiana path.
async function loadApiMetadata() {
    try {
        const response = await fetch(API_ENDPOINT, { method: 'GET' });
        if (response.ok) {
            const data = await response.json();
            if (data.host_metadata && typeof data.host_metadata === 'object') {
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

function initEventListeners() {
    // Input Handling
    elements.fileUpload.addEventListener('change', handleFileUpload);
    elements.sequenceInput.addEventListener('input', debounce(handleSequenceChange, 300));
    elements.clearBtn.addEventListener('click', clearAll);

    // Objective Change
    elements.objectiveRadios.forEach(radio => {
        radio.addEventListener('change', (e) => {
            state.objective = e.target.value;
            updateDesignBriefSummary();
        });
    });
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
        updateDesignBriefSummary();
    });

    elements.kozakToggle.addEventListener('change', (e) => {
        state.kozak = e.target.checked;
        updateDesignBriefSummary();
    });
    elements.dinucToggle.addEventListener('change', (e) => {
        state.dinuc = e.target.checked;
        updateDesignBriefSummary();
    });
    elements.customRestrictionSites.addEventListener('input', () => {
        state.customRestrictionSites = [];
        updateDesignBriefSummary();
    });
    elements.typeIisEnzymes.forEach(input => input.addEventListener('change', updateDesignBriefSummary));
    elements.saveReviewerDisposition.addEventListener('click', saveReviewerDisposition);
    elements.clearHistory.addEventListener('click', clearHistory);

    // Action
    elements.optimizeBtn.addEventListener('click', runOptimization);

    // Results Actions
    elements.downloadFasta.addEventListener('click', () => downloadFile('fasta'));
    elements.downloadGenbank.addEventListener('click', () => downloadFile('genbank'));
    elements.copyBtn.addEventListener('click', copyToClipboard);
    elements.copyConstructId.addEventListener('click', copyConstructId);
    elements.submitValidationBtn.addEventListener('click', submitValidation);
    elements.copyJsonBtn.addEventListener('click', copyJson);
    elements.toggleDetails.addEventListener('click', toggleDetailsPanel);
    elements.themeToggle.addEventListener('click', toggleTheme);
    elements.changelogBtn.addEventListener('click', toggleChangelog);
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
        feasibility_best: 'recommended feasibility design',
        high_cai: 'CAI-focused comparison',
        gc_target: 'GC-focused comparison',
        assembly_friendly: 'assembly-oriented comparison'
    };
    const requirements = [];
    if (elements.useTemplateCheck?.checked) requirements.push('MoClo');
    const enzymes = Array.from(elements.typeIisEnzymes || [])
        .filter(input => input.checked)
        .map(input => input.value);
    if (enzymes.length) requirements.push(`avoid ${enzymes.join('/')}`);
    if (elements.customRestrictionSites?.value.trim()) requirements.push('custom restriction sites');
    if (elements.kozakToggle?.checked) requirements.push('Kozak handling');
    if (elements.dinucToggle?.checked) requirements.push('TpA reduction');
    const requirementText = requirements.length ? ` · ${requirements.join(' · ')}` : '';
    const method = state.host === 'by2' && selectedObjective === 'feasibility_best'
        ? 'stable profile design'
        : methodLabels[selectedObjective] || 'deterministic design';
    elements.appliedPolicySummary.textContent = `${host} · ${method}${requirementText}`;
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

async function runOptimization() {
    // Flush the debounced input handler so an immediate click after paste/type
    // optimizes the current textarea value instead of stale state.
    handleSequenceChange({ target: elements.sequenceInput });

    if (!state.sequence || state.sequence.length < 3) {
        showToast('Please enter a valid DNA sequence', 'error');
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
        // Prepare Request
        const payload = {
            sequence: state.sequence,
            host: state.host,
            mode: state.engineMode,
            save_db: state.saveDb,
            use_template: state.useTemplate,
            kozak: state.kozak,
            dinuc: state.dinuc,
            return_candidates: true
        };
        if (state.engineMode !== 'profile') {
            payload.profile = 'balanced';
        } else if (state.objective === 'feasibility_best' && state.host === 'nbenthamiana') {
            payload.objective = 'feasibility_best';
            payload.host_profile = state.host;
            // Host-aware default (v3.3.0) — previously hardcoded to
            // the legacy 55-65 band, which silently overrode the server's
            // resolve_host_gc_range() default for every feasibility_best run.
            payload.constraints = getGcRange(state.host);
        } else {
            payload.profile = state.objective === 'feasibility_best' ? 'balanced' : state.objective;
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

        if (ENABLE_MOCK) {
            showToast('Running mock optimization for demonstration...', 'info');
            await new Promise(r => setTimeout(r, 1500));
            state.results = getMockResult();
            renderResults();
            showToast('Showing simulated results', 'success');
        } else {
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

function renderResults() {
    const res = state.results;
    if (!res) return;
    const primary = getPrimaryResult(res);

    elements.emptyState.classList.add('hidden');
    elements.resultsContainer.classList.remove('hidden');
    if (elements.resultContextSummary) {
        const host = formatHostProfile(getResultHostProfile(res));
        elements.resultContextSummary.textContent = `${host} · Review the computational checks before synthesis or experimental use.`;
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
        renderResultsReport(res, primary, gcTarget);
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

    const recommendedId = res.recommended_candidate ? res.recommended_candidate.id : res.candidates[0].id;
    elements.candidateComparisonBody.innerHTML = res.candidates.map(candidate => {
        const isRecommended = candidate.id === recommendedId;
        const status = candidate.automated_decision || (candidate.validator_status === 'pass' ? 'PASS' : 'REVIEW');
        const aaStatus = candidate.validator_status === 'pass' ? '100%' : 'Review';
        const description = CANDIDATE_DESCRIPTIONS[candidate.id];
        const nameCell = description
            ? `${escapeHtml(candidate.label || candidate.id)}<span class="tooltip-container"><span class="info-icon">i</span><span class="tooltip-text">${escapeHtml(description)}</span></span>`
            : escapeHtml(candidate.label || candidate.id);
        return `
            <tr class="${isRecommended ? 'bg-emerald-50 dark:bg-emerald-900/20' : 'bg-white dark:bg-slate-900'}">
                <td class="px-3 py-2 font-bold text-slate-800 dark:text-slate-100">${nameCell}${isRecommended ? ' ★' : ''}</td>
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
    elements.mfeWarningBanner.textContent = `RNA secondary-structure/MFE analysis was not performed (${mfe.reason}).`;
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
    const criteriaRows = checks.filter(check => check.id !== 'mfe');
    const requiredFailures = Number.isInteger(summary.required_failure_count)
        ? summary.required_failure_count
        : (criteriaRows.length ? criteriaRows.filter(check => check.status === 'FAIL').length : null);
    const preferredWarnings = Number.isInteger(summary.preferred_warning_count)
        ? summary.preferred_warning_count
        : (criteriaRows.length ? criteriaRows.filter(check => check.status === 'WARNING').length : null);

    const reviewPlan = buildResearcherReviewPlan(checks, decision, res.seed ?? null, settingMismatches);
    return {
        report_schema_version: '1.0',
        identity: {
            result_id: res.result_identifier || res.construct_id || null,
            construct_id: res.construct_id || null,
            result_created_at: res.created_at || null,
            report_generated_at: new Date().toISOString(),
        },
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
        artifacts: { optimized_sequence: outputSequence },
    };
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

function renderResultsReport(res, primary, gcTarget) {
    if (!elements.resultsReport || !elements.resultsReportBody) return;
    const model = buildResultsReportModel(res, primary, gcTarget);
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
    const { artifacts: _artifacts, ...evidence } = model;
    return evidence;
}

function standaloneReportHtml(model) {
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
    const provenanceRows = [
        ['Product version', model.provenance.product_version], ['Engine', model.context.engine], ['Objective', model.context.objective],
        ['Host / profile', `${model.context.host_profile} / ${model.context.profile}`], ['Codon reference', model.provenance.codon_reference_id],
        ['Reference policy', model.provenance.reference_policy_version], ['GC reference band', model.provenance.gc_reference_band],
        ['Seed', model.context.seed === null ? 'Not specified' : model.context.seed], ['Result created (UTC)', model.identity.result_created_at],
        ['Input SHA-256', model.provenance.input_sequence_hash], ['Output SHA-256', model.provenance.output_cds_hash], ['Parameter SHA-256', model.provenance.parameter_hash],
    ].map(([label, value]) => `<dt>${escapeHtml(label)}</dt><dd>${escapeHtml(reportValue(value))}</dd>`).join('');
    const sequence = model.artifacts.optimized_sequence.match(/.{1,60}/g)?.join('\n') || model.artifacts.optimized_sequence;
    return `<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>FactorForge Researcher Decision Report — ${escapeHtml(reportValue(model.identity.result_id))}</title><style>
body{margin:0;background:#f1f5f9;color:#0f172a;font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}.page{max-width:920px;margin:auto;padding:42px 24px 72px}h1{font-size:30px;margin:4px 0}h2{font-size:14px;text-transform:uppercase;letter-spacing:.08em;color:#475569;margin:28px 0 10px}h3{margin:0}.eyebrow{color:#047857;font-weight:800;text-transform:uppercase;letter-spacing:.08em}.muted,small{color:#64748b}.outcome{background:#fff;border:2px solid;border-radius:16px;padding:20px}.outcome strong{display:block;font-size:28px}.counts{font-weight:700}.callout{padding:14px;border-radius:12px;background:#fff7ed;border:1px solid #fdba74}.priority{display:grid;grid-template-columns:32px 1fr;gap:10px;background:#fff;border:1px solid #dbe2ea;border-left:4px solid;border-radius:12px;padding:14px;margin:8px 0}.priority-number{width:26px;height:26px;border-radius:50%;background:#e2e8f0;display:grid;place-items:center;font-weight:800}.scroll{overflow-x:auto}table{width:100%;border-collapse:collapse;background:#fff}th,td{text-align:left;vertical-align:top;padding:9px;border-bottom:1px solid #e2e8f0}.mismatch{background:#fff1f2;color:#9f1239}dl{display:grid;grid-template-columns:minmax(150px,220px) 1fr;gap:6px 14px}dt{font-weight:700}dd{margin:0;font-family:monospace;overflow-wrap:anywhere}pre{background:#0f172a;color:#a7f3d0;padding:16px;border-radius:12px;overflow:auto;font:12px/1.6 monospace}ul{padding-left:20px}@media(prefers-color-scheme:dark){body{background:#0b0f19;color:#f1f5f9}h2{color:#94a3b8}.eyebrow{color:#34d399}.muted,small{color:#94a3b8}.outcome{background:#1e293b;border-color:#334155}.callout{background:#451a03;border-color:#9a3412;color:#fed7aa}.priority{background:#1e293b;border-color:#334155}.priority-number{background:#334155;color:#f8fafc}table{background:#1e293b}th,td{border-bottom:1px solid #334155}.mismatch{background:#4c0519;color:#fda4af}pre{background:#020617;color:#6ee7b7;border:1px solid #1e293b}}@media(max-width:640px){.page{padding:24px 14px}dl{grid-template-columns:1fr}dd{margin-bottom:7px}}@media print{body{background:#fff!important;color:#0f172a!important}.page{padding:0}.outcome,.priority,table{background:#fff!important;break-inside:avoid}.callout{border-color:#999}}
</style></head><body><main class="page"><p class="eyebrow">FactorForge CDS Design Review</p><h1>Researcher Decision Report</h1><p class="muted">Result ${escapeHtml(reportValue(model.identity.result_id))} · created ${escapeHtml(reportValue(model.identity.result_created_at))} · report generated ${escapeHtml(model.identity.report_generated_at)}</p><div class="callout"><b>Sequence-data notice:</b> This HTML contains the optimized DNA sequence. Handle and share it according to your sequence-data policy.</div><section><h2>Decision brief</h2><div class="outcome" style="border-color:${decisionColor}"><strong style="color:${decisionColor}">${escapeHtml(model.disposition.automated_decision.replaceAll('_', ' '))}</strong><p>${escapeHtml(model.interpretation.headline)}</p><p class="counts">${escapeHtml(`${reportValue(model.disposition.required_failure_count)} required fail · ${reportValue(model.disposition.preferred_warning_count)} warning · ${reportValue(model.disposition.unavailable_check_count)} unavailable`)}</p><small>${escapeHtml(reportValue(model.disposition.explanation))}</small></div></section><section><h2>Review priorities and next actions</h2>${priorities}</section><section><h2>Requested vs applied settings</h2><div class="scroll"><table><thead><tr><th>Setting</th><th>Requested</th><th>Applied / recorded</th></tr></thead><tbody>${settings}</tbody></table></div></section><section><h2>All computational checks</h2><div class="scroll"><table><thead><tr><th>Check</th><th>Observed</th><th>Policy</th><th>Status</th></tr></thead><tbody>${checks}</tbody></table></div></section>${candidates}<section><h2>Sequence and process</h2><p>Input type: <b>${escapeHtml(model.context.input_type)}</b> · output length: <b>${escapeHtml(reportValue(model.sequence_summary.output_length_nt, ' nt'))}</b> · nucleotide comparison: <b>${escapeHtml(model.sequence_summary.comparison_available ? reportValue(model.sequence_summary.nucleotide_changes) : 'Not recorded')}</b></p><p>Domestication: ${escapeHtml(model.process.domestication_attempted ? 'Attempted' : 'Not attempted')} · MFE: ${escapeHtml(model.metrics.mfe_status === 'computed' ? 'Computed' : `Not computed (${model.metrics.mfe_status_reason})`)}</p></section><section><h2>Downstream handoff checklist</h2><ul>${model.interpretation.next_steps.map(item => `<li>${escapeHtml(item)}</li>`).join('')}</ul></section><section><h2>Reproducibility and provenance</h2><dl>${provenanceRows}</dl></section><section><h2>Optimized sequence (DNA)</h2><pre>${escapeHtml(sequence)}</pre></section><section><h2>Interpretation and limitations</h2><p>${escapeHtml(model.interpretation.scope)}</p><ul>${model.interpretation.limitations.map(item => `<li>${escapeHtml(item)}</li>`).join('')}</ul></section></main></body></html>`;
}

function downloadResultsReportHtml(model) {
    trackEvent('report_download', { format: 'html' });
    downloadTextArtifact(standaloneReportHtml(model), 'text/html;charset=utf-8', `factorforge_design_review_${reportFileStem(model)}.html`);
    showToast('Design review report downloaded', 'success');
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
        automated_decision: candidate.automated_decision || null,
        required_failure_count: candidate.required_failure_count ?? null,
        preferred_warning_count: candidate.preferred_warning_count ?? null,
    })) : [];
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
        profile: state.objective,
        host: getResultHostProfile(result),
        cai: primary.metrics.cai,
        gc: calculateGC(primary.optimized_sequence),
        sequence: primary.optimized_sequence,
        resultSnapshot: historyResultSnapshot(result, primary),
    };

    state.history.unshift(item);
    if (state.history.length > 10) state.history.pop();
    localStorage.setItem('factorforge_history', JSON.stringify({ schemaVersion: HISTORY_SCHEMA_VERSION, items: state.history }));
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

    if (format === 'fasta') {
        content = `>FactorForge_Optimized | Objective: ${state.objective} | CAI: ${primary.metrics.cai}\n${seq}`;
        fileName = `optimized_sequence_${Date.now()}.fasta`;
    } else {
        // Basic GenBank template
        content = `LOCUS       Exported                ${seq.length} bp    DNA     linear   \n`;
        const hostLabel = HOST_LABELS[state.host] || 'N. benthamiana';
        content += `DEFINITION  FactorForge Optimized Sequence for ${hostLabel}\n`;
        content += `FEATURES             Location/Qualifiers\n`;
        content += `     CDS             1..${seq.length}\n`;
        content += `                     /label="Optimized_CDS"\n`;
        content += `                     /note="Objective: ${state.objective}"\n`;
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
        const version = state.results.engine_versions?.product || '3.4.5';
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
