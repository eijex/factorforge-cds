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
const HISTORY_SCHEMA_VERSION = 2;

function loadVersionedHistory() {
    try {
        const raw = JSON.parse(localStorage.getItem('factorforge_history') || '[]');
        if (Array.isArray(raw)) return raw.map(item => ({ ...item, schemaVersion: HISTORY_SCHEMA_VERSION }));
        if (raw && Array.isArray(raw.items)) return raw.items;
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
    hostSelect: document.getElementById('hostSelect'),
    saveDbToggle: document.getElementById('saveDbToggle'),
    saveDbStatus: document.getElementById('saveDbStatus'),
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
            const dbAvailable = Boolean(apiCapabilities.db_save?.available);
            elements.saveDbToggle.disabled = !dbAvailable;
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
    });
    elements.saveDbToggle.addEventListener('change', (e) => {
        state.saveDb = e.target.checked;
    });
    elements.alignmentPrev.addEventListener('click', () => changeAlignmentPage(-1));
    elements.alignmentNext.addEventListener('click', () => changeAlignmentPage(1));

    elements.useTemplateCheck.addEventListener('change', (e) => {
        state.useTemplate = e.target.checked;
    });

    elements.kozakToggle.addEventListener('change', (e) => state.kozak = e.target.checked);
    elements.dinucToggle.addEventListener('change', (e) => state.dinuc = e.target.checked);
    elements.customRestrictionSites.addEventListener('input', () => {
        state.customRestrictionSites = [];
    });
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
    if (!seq || seq.length < 10) return;
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
    const optSeq = primary.optimized_sequence;
    const calculatedGC = calculateGC(optSeq);
    const oGC = isProteinInput ? null : calculateGC(origSeq);

    elements.origGC.textContent = isProteinInput ? 'N/A' : `${oGC}%`;
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
    return {
        status: metrics.mfe_status || 'not_computed',
        reason: metrics.mfe_status_reason || 'status unavailable'
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

function reportStatCard({ label, value, sub, tone = 'neutral' }) {
    const toneClasses = {
        good: 'bg-emerald-50 dark:bg-emerald-900/20 text-emerald-700 dark:text-emerald-300',
        warn: 'bg-amber-50 dark:bg-amber-900/20 text-amber-700 dark:text-amber-300',
        bad: 'bg-rose-50 dark:bg-rose-900/20 text-rose-700 dark:text-rose-300',
        neutral: 'bg-slate-100 dark:bg-slate-800 text-slate-700 dark:text-slate-200',
    }[tone] || 'bg-slate-100 dark:bg-slate-800 text-slate-700 dark:text-slate-200';

    return `
        <div class="min-w-0 rounded-xl p-2.5 ${toneClasses}">
            <p class="text-[9px] font-extrabold uppercase tracking-widest opacity-70 truncate">${escapeHtml(label)}</p>
            <p class="mt-1 text-xs font-black leading-tight break-words">${escapeHtml(value)}</p>
            ${sub ? `<p class="mt-0.5 text-[10px] font-medium opacity-80 break-words">${escapeHtml(sub)}</p>` : ''}
        </div>
    `;
}

function computeResultsReportData(res, primary, gcTarget) {
    const host = getResultHostProfile(res);
    const profile = res.profile || state.objective;
    const seedText = res.seed == null ? 'seed not specified' : `seed=${res.seed}`;

    const cai = Number(primary.metrics.cai || 0);
    const gc = Number(primary.metrics.gc_percent || 0);
    const gcInRange = gc >= gcTarget.min && gc <= gcTarget.max;

    const custom = res.custom_restriction_sites;
    const removed = Array.isArray(custom?.removed) ? custom.removed : [];
    const unresolved = Array.isArray(custom?.unresolved) ? custom.unresolved : [];
    const attempted = Boolean(custom);
    const requestedNames = Array.isArray(custom?.requested)
        ? custom.requested.map(site => site.name).filter(name => Object.hasOwn(TYPE_IIS_PRESETS, name))
        : [];
    const selected = requestedNames.length > 0 ? requestedNames : state.selectedTypeIisEnzymes;
    const typeIisFail = unresolved.length > 0;

    const mfe = getMfeStatus(res);
    const mfeComputed = mfe.status === 'computed';

    return {
        generatedAt: new Date().toISOString(),
        host, profile, seedText, seed: res.seed ?? null,
        cai, gc, gcTarget, gcInRange,
        typeIis: { selected, fail: typeIisFail, checked: selected.length > 0 },
        domestication: { attempted, removed, unresolved },
        mfe: { computed: mfeComputed, reason: mfe.reason },
        candidates: Array.isArray(res.candidates) ? res.candidates : [],
    };
}

function reportCardsFromData(data) {
    return [
        reportStatCard({ label: 'Host / Profile', value: `${data.host} · ${data.profile}`, sub: data.seedText, tone: 'neutral' }),
        reportStatCard({
            label: 'CAI',
            value: data.cai.toFixed(3),
            sub: data.cai >= 0.8 ? 'meets 0.800 minimum' : 'below 0.800 minimum',
            tone: data.cai >= 0.8 ? 'good' : 'warn',
        }),
        reportStatCard({
            label: 'GC content',
            value: `${data.gc.toFixed(1)}%`,
            sub: `target ${data.gcTarget.min.toFixed(1)}–${data.gcTarget.max.toFixed(1)}%`,
            tone: data.gcInRange ? 'good' : 'warn',
        }),
        reportStatCard({
            label: 'Type IIS',
            value: data.typeIis.checked ? (data.typeIis.fail ? 'FAIL' : 'PASS') : 'Not checked',
            sub: data.typeIis.checked ? data.typeIis.selected.join(', ') : 'no preset enzymes selected',
            tone: !data.typeIis.checked ? 'neutral' : (data.typeIis.fail ? 'bad' : 'good'),
        }),
        reportStatCard({
            label: 'Domestication',
            value: data.domestication.attempted ? 'Attempted' : 'Not attempted',
            sub: data.domestication.attempted
                ? `${data.domestication.removed.length} removed · ${data.domestication.unresolved.length} unresolved`
                : 'no enzymes selected to fix',
            tone: !data.domestication.attempted ? 'neutral' : (data.domestication.unresolved.length > 0 ? 'warn' : 'good'),
        }),
        reportStatCard({
            label: 'MFE',
            value: data.mfe.computed ? 'Computed' : 'Not computed',
            sub: data.mfe.computed ? '' : data.mfe.reason,
            tone: data.mfe.computed ? 'good' : 'warn',
        }),
    ];
}

function renderResultsReport(res, primary, gcTarget) {
    if (!elements.resultsReport || !elements.resultsReportBody) return;
    const data = computeResultsReportData(res, primary, gcTarget);
    const cards = reportCardsFromData(data);

    const comparisonRows = data.candidates.length > 1
        ? `
            <div class="overflow-x-auto pt-1">
                <p class="text-[10px] font-extrabold uppercase tracking-widest text-slate-500 dark:text-slate-400 mb-2">Candidate comparison</p>
                <table class="w-full text-left text-xs">
                    <thead>
                        <tr class="text-slate-500 dark:text-slate-400">
                            <th class="py-1 pr-3 font-bold">Profile</th>
                            <th class="py-1 pr-3 font-bold">CAI</th>
                            <th class="py-1 font-bold">GC%</th>
                        </tr>
                    </thead>
                    <tbody class="divide-y divide-slate-100 dark:divide-slate-800">
                        ${data.candidates.map(candidate => `
                            <tr>
                                <td class="py-1.5 pr-3 font-semibold text-slate-700 dark:text-slate-200">${escapeHtml(candidate.label || candidate.id)}</td>
                                <td class="py-1.5 pr-3 font-mono">${Number(candidate.cai || 0).toFixed(3)}</td>
                                <td class="py-1.5 font-mono">${Number(candidate.gc_percent || 0).toFixed(1)}</td>
                            </tr>
                        `).join('')}
                    </tbody>
                </table>
            </div>
        `
        : '';

    elements.resultsReportBody.innerHTML = `
        <div class="grid grid-cols-2 gap-2">${cards.join('')}</div>
        ${comparisonRows}
        <button type="button" id="downloadResultsReportBtn" class="mt-3 w-full flex items-center justify-center gap-2 rounded-xl bg-slate-800 hover:bg-slate-700 dark:bg-slate-700 dark:hover:bg-slate-600 text-white text-xs font-bold py-2.5 transition-colors">
            <span>📄</span> Download Report (HTML)
        </button>
    `;
    elements.resultsReport.classList.remove('hidden');

    const downloadBtn = document.getElementById('downloadResultsReportBtn');
    if (downloadBtn) {
        downloadBtn.addEventListener('click', () => downloadResultsReportHtml(res, primary, gcTarget));
    }
}

function downloadResultsReportHtml(res, primary, gcTarget) {
    trackEvent('report_download', { format: 'html' });
    const data = computeResultsReportData(res, primary, gcTarget);
    const toneHex = { good: '#059669', warn: '#d97706', bad: '#e11d48', neutral: '#475569' };
    const cardTone = (t) => toneHex[t] || toneHex.neutral;

    const cardDefs = [
        { label: 'Host / Profile', value: `${data.host} · ${data.profile}`, sub: data.seedText, tone: 'neutral' },
        { label: 'CAI', value: data.cai.toFixed(3), sub: data.cai >= 0.8 ? 'meets 0.800 minimum' : 'below 0.800 minimum', tone: data.cai >= 0.8 ? 'good' : 'warn' },
        { label: 'GC content', value: `${data.gc.toFixed(1)}%`, sub: `target ${data.gcTarget.min.toFixed(1)}–${data.gcTarget.max.toFixed(1)}%`, tone: data.gcInRange ? 'good' : 'warn' },
        { label: 'Type IIS', value: data.typeIis.checked ? (data.typeIis.fail ? 'FAIL' : 'PASS') : 'Not checked', sub: data.typeIis.checked ? data.typeIis.selected.join(', ') : 'no preset enzymes selected', tone: !data.typeIis.checked ? 'neutral' : (data.typeIis.fail ? 'bad' : 'good') },
        { label: 'Domestication', value: data.domestication.attempted ? 'Attempted' : 'Not attempted', sub: data.domestication.attempted ? `${data.domestication.removed.length} removed · ${data.domestication.unresolved.length} unresolved` : 'no enzymes selected to fix', tone: !data.domestication.attempted ? 'neutral' : (data.domestication.unresolved.length > 0 ? 'warn' : 'good') },
        { label: 'MFE', value: data.mfe.computed ? 'Computed' : 'Not computed', sub: data.mfe.computed ? '' : data.mfe.reason, tone: data.mfe.computed ? 'good' : 'warn' },
    ];

    const cardsHtml = cardDefs.map(c => `
        <div style="border-radius:14px;padding:18px;background:#f8fafc;border:1px solid #e2e8f0;border-left:4px solid ${cardTone(c.tone)};">
            <p style="margin:0;font-size:11px;font-weight:800;letter-spacing:.08em;text-transform:uppercase;color:#64748b;">${escapeHtml(c.label)}</p>
            <p style="margin:6px 0 0;font-size:20px;font-weight:800;color:${cardTone(c.tone)};">${escapeHtml(c.value)}</p>
            ${c.sub ? `<p style="margin:4px 0 0;font-size:12px;color:#475569;">${escapeHtml(c.sub)}</p>` : ''}
        </div>
    `).join('');

    const comparisonHtml = data.candidates.length > 1 ? `
        <h2 style="font-size:14px;text-transform:uppercase;letter-spacing:.08em;color:#64748b;margin:32px 0 12px;">Candidate comparison</h2>
        <table style="width:100%;border-collapse:collapse;font-size:13px;">
            <thead><tr style="text-align:left;color:#64748b;">
                <th style="padding:6px 12px 6px 0;border-bottom:1px solid #e2e8f0;">Profile</th>
                <th style="padding:6px 12px 6px 0;border-bottom:1px solid #e2e8f0;">CAI</th>
                <th style="padding:6px 0;border-bottom:1px solid #e2e8f0;">GC%</th>
            </tr></thead>
            <tbody>
                ${data.candidates.map(c => `
                    <tr>
                        <td style="padding:8px 12px 8px 0;border-bottom:1px solid #f1f5f9;font-weight:600;">${escapeHtml(c.label || c.id)}</td>
                        <td style="padding:8px 12px 8px 0;border-bottom:1px solid #f1f5f9;font-family:monospace;">${Number(c.cai || 0).toFixed(3)}</td>
                        <td style="padding:8px 0;border-bottom:1px solid #f1f5f9;font-family:monospace;">${Number(c.gc_percent || 0).toFixed(1)}</td>
                    </tr>
                `).join('')}
            </tbody>
        </table>
    ` : '';

    const seq = primary.optimized_sequence || '';
    const seqWrapped = seq.match(/.{1,60}/g)?.join('\n') || seq;

    const html = `<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>FactorForge Results Report — ${escapeHtml(data.host)} / ${escapeHtml(data.profile)}</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
  body { margin:0; background:#f1f5f9; color:#0f172a; font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; }
  .page { max-width:820px; margin:0 auto; padding:48px 24px 80px; }
  .grid { display:grid; grid-template-columns:repeat(3, 1fr); gap:14px; }
  @media (max-width:640px) { .grid { grid-template-columns:repeat(2, 1fr); } }
  pre { background:#0f172a; color:#a7f3d0; padding:16px; border-radius:12px; overflow-x:auto; font-size:12px; line-height:1.6; }
</style>
</head>
<body>
<div class="page">
  <p style="font-size:12px;font-weight:700;letter-spacing:.08em;text-transform:uppercase;color:#059669;margin:0;">FactorForge CDS Design Review</p>
  <h1 style="font-size:26px;margin:6px 0 4px;">Results Report</h1>
  <p style="color:#64748b;font-size:13px;margin:0 0 28px;">Generated ${escapeHtml(new Date(data.generatedAt).toLocaleString())}</p>
  <div class="grid">${cardsHtml}</div>
  ${comparisonHtml}
  <h2 style="font-size:14px;text-transform:uppercase;letter-spacing:.08em;color:#64748b;margin:32px 0 12px;">Optimized sequence (DNA)</h2>
  <pre>${escapeHtml(seqWrapped)}</pre>
  <p style="margin-top:32px;font-size:11px;color:#94a3b8;">This is an in-silico CDS design candidate and pre-synthesis review artifact. Review and wet-lab testing are required before relying on any design in experiments.</p>
</div>
</body>
</html>`;

    const blob = new Blob([html], { type: 'text/html' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `factorforge_results_report_${Date.now()}.html`;
    a.click();
    window.URL.revokeObjectURL(url);
    showToast('Report downloaded', 'success');
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

function addToHistory(input, result) {
    const primary = getPrimaryResult(result);
    const item = {
        id: Date.now(),
        timestamp: new Date().toLocaleString(),
        inputLen: input.length,
        profile: state.objective,
        host: getResultHostProfile(result),
        cai: primary.metrics.cai,
        gc: calculateGC(primary.optimized_sequence),
        sequence: primary.optimized_sequence,
        inputSequence: input
        ,acceptanceCriteria: result.acceptance_criteria_snapshot || null
        ,automatedDecision: result.automated_decision || null
        ,reviewerDisposition: result.reviewer_disposition || null
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
                <p class="text-[10px] font-bold text-slate-700">${item.inputLen}bp → ${item.profile}</p>
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

    elements.sequenceInput.value = item.inputSequence;
    handleSequenceChange({ target: { value: item.inputSequence } });
    state.results = {
        optimized_sequence: item.sequence,
        metrics: { cai: item.cai, gc_percent: item.gc, polya_signals: 0, length: item.sequence.length },
        profile: item.profile,
        host_profile: item.host || 'nbenthamiana',
        acceptance_criteria_snapshot: item.acceptanceCriteria || null,
        automated_decision: item.automatedDecision || null,
        reviewer_disposition: item.reviewerDisposition || null
    };
    renderResults();
    showToast('History item loaded', 'success');
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
