/**
 * FactorForge Application Logic
 * Vanilla JS for the FactorForge CDS design interface.
 */

const API_ENDPOINT = '/api/optimize';
const ENABLE_MOCK = window.FACTORFORGE_ENABLE_MOCK === true;
const HOST_LABELS = {
    nbenthamiana: 'N. benthamiana',
    by2: 'Tobacco BY-2',
    ntabacum: 'Tobacco BY-2'
};

// State Management
const state = {
    sequence: '',
    objective: 'feasibility_best',
    host: 'nbenthamiana',
    useTemplate: false,
    kozak: false,
    dinuc: false,
    customRestrictionSites: [],
    results: null,
    isOptimizing: false,
    history: JSON.parse(localStorage.getItem('factorforge_history') || '[]')
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
    hostRadios: document.getElementsByName('host'),
    feasibilityBestOption: document.getElementById('feasibilityBestOption'),
    feasibilityBestCard: document.getElementById('feasibilityBestCard'),
    feasibilityBestHostBadge: document.getElementById('feasibilityBestHostBadge'),
    useTemplateCheck: document.getElementById('useTemplate'),
    kozakToggle: document.getElementById('toggleKozak'),
    dinucToggle: document.getElementById('toggleDinuc'),
    customRestrictionSites: document.getElementById('customRestrictionSites'),
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
    gcChart: document.getElementById('gcChart'),
    historyList: document.getElementById('historyList'),
    clearHistory: document.getElementById('clearHistory'),
    changelogBtn: document.getElementById('changelogBtn'),
    changelogModal: document.getElementById('changelogModal'),
    closeModal: document.getElementById('closeModal'),
    modalOverlay: document.getElementById('modalOverlay'),
    inputTypeBadge: document.getElementById('inputTypeBadge'),
    toastContainer: document.getElementById('toastContainer'),
    logoIcon: document.getElementById('logoIcon'),
    logoTitle: document.getElementById('logoTitle')
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
function gcBucket(gc) {
    if (gc < 40) return '<40';
    if (gc < 55) return '40-55';
    if (gc < 65) return '55-65';
    return '>65';
}
function trackEvent(name, data) {
    try { window.va?.('event', { name, data }); } catch (_) {}
}

// Initialization
document.addEventListener('DOMContentLoaded', () => {
    initTheme();
    applyStaticLabelPatches();
    initEventListeners();
    updateHostUI();
    renderHistory();
    console.log('FactorForge v3.2.0 Engaged');
});

function initEventListeners() {
    // Input Handling
    elements.fileUpload.addEventListener('change', handleFileUpload);
    elements.sequenceInput.addEventListener('input', debounce(handleSequenceChange, 300));
    elements.clearBtn.addEventListener('click', clearAll);

    // Objective Change
    elements.objectiveRadios.forEach(radio => {
        radio.addEventListener('change', (e) => {
            state.objective = e.target.value;
            updateHostUI();
        });
    });

    elements.hostRadios.forEach(radio => {
        radio.addEventListener('change', (e) => {
            state.host = e.target.value;
            updateHostUI();
        });
    });

    elements.useTemplateCheck.addEventListener('change', (e) => {
        state.useTemplate = e.target.checked;
    });

    elements.kozakToggle.addEventListener('change', (e) => state.kozak = e.target.checked);
    elements.dinucToggle.addEventListener('change', (e) => state.dinuc = e.target.checked);
    elements.customRestrictionSites.addEventListener('input', () => {
        state.customRestrictionSites = [];
    });
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
}

function reloadPage() {
    window.location.reload();
}

function applyStaticLabelPatches() {
    if (!elements.submitValidationBtn) return;
    const label = elements.submitValidationBtn.children[1];
    if (label) label.textContent = 'Share Wet-Lab Result';
}

function updateHostUI() {
    const selectedHost = Array.from(elements.hostRadios).find(radio => radio.checked);
    state.host = selectedHost?.value || state.host || 'nbenthamiana';

    const isBy2 = state.host === 'by2';
    const objectiveRadios = Array.from(elements.objectiveRadios);
    const feasibilityRadio = objectiveRadios.find(radio => radio.value === 'feasibility_best');

    if (feasibilityRadio) {
        feasibilityRadio.disabled = isBy2;
        if (isBy2 && feasibilityRadio.checked) {
            const fallback = objectiveRadios.find(radio => !radio.disabled && radio.value === 'high_cai')
                || objectiveRadios.find(radio => !radio.disabled);
            if (fallback) {
                fallback.checked = true;
                state.objective = fallback.value;
            }
        }
    }

    if (elements.feasibilityBestOption) {
        elements.feasibilityBestOption.classList.toggle('cursor-not-allowed', isBy2);
    }
    if (elements.feasibilityBestCard) {
        elements.feasibilityBestCard.classList.toggle('opacity-60', isBy2);
        elements.feasibilityBestCard.classList.toggle('cursor-not-allowed', isBy2);
        elements.feasibilityBestCard.classList.toggle('cursor-pointer', !isBy2);
        elements.feasibilityBestCard.classList.toggle('hover:border-emerald-400', !isBy2);
        elements.feasibilityBestCard.classList.toggle('hover:shadow-md', !isBy2);
    }
    if (elements.feasibilityBestHostBadge) {
        elements.feasibilityBestHostBadge.classList.toggle('hidden', !isBy2);
    }
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
    if (state.results) renderGCGraph(state.results.optimized_sequence);
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

    // Clean sequence (remove FASTA headers, whitespace)
    const lines = rawValue.split('\n');
    let sequenceOnly = lines[0].startsWith('>') ? lines.slice(1).join('') : rawValue;
    sequenceOnly = sequenceOnly.replace(/[^A-Za-z*]/g, ''); // Include * for stop codons

    // Detection Regex
    const dnaRegex = /^[ACGTacgt]+$/;
    const proteinRegex = /^[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy*]+$/;

    const isDNA = dnaRegex.test(sequenceOnly);
    const isProtein = proteinRegex.test(sequenceOnly);

    // Update Badge & Warning
    const dnaInputWarning = document.getElementById('dnaInputWarning');
    if (sequenceOnly.length > 0) {
        elements.inputTypeBadge.classList.remove('hidden');
        if (isDNA) {
            elements.inputTypeBadge.textContent = 'DNA';
            elements.inputTypeBadge.className = 'ml-2 px-1.5 py-0.5 rounded bg-blue-100 text-blue-700 font-bold text-[9px] uppercase';
            elements.validationWarning.classList.add('hidden');
            if (dnaInputWarning) dnaInputWarning.classList.remove('hidden');
        } else if (isProtein) {
            elements.inputTypeBadge.textContent = 'Protein';
            elements.inputTypeBadge.className = 'ml-2 px-1.5 py-0.5 rounded bg-emerald-100 text-emerald-700 font-bold text-[9px] uppercase';
            elements.validationWarning.classList.add('hidden');
            if (dnaInputWarning) dnaInputWarning.classList.add('hidden');
        } else {
            elements.inputTypeBadge.textContent = 'Mixed/Invalid';
            elements.inputTypeBadge.className = 'ml-2 px-1.5 py-0.5 rounded bg-rose-100 text-rose-700 font-bold text-[9px] uppercase';
            elements.validationWarning.classList.remove('hidden');
            elements.validationWarning.innerHTML = '<span class="mr-2">⚠️</span> Contains characters not valid for DNA or Protein.';
            if (dnaInputWarning) dnaInputWarning.classList.add('hidden');
        }
    } else {
        elements.inputTypeBadge.classList.add('hidden');
        elements.validationWarning.classList.add('hidden');
        if (dnaInputWarning) dnaInputWarning.classList.add('hidden');
    }

    state.sequence = sequenceOnly.toUpperCase();

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
        elements.inputLenBadge.style.opacity = "0.3";
        elements.inputGCBadge.textContent = "GC: 0%";
        elements.inputGCBadge.style.opacity = "0.3";
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
            use_template: state.useTemplate,
            kozak: state.kozak,
            dinuc: state.dinuc,
            return_candidates: true
        };
        if (state.objective === 'feasibility_best') {
            payload.objective = 'feasibility_best';
            payload.host_profile = state.host;
            payload.constraints = { gc_min: 55.0, gc_max: 65.0 };
        } else {
            payload.profile = state.objective;
        }
        const customRestrictionSites = parseCustomRestrictionSites(elements.customRestrictionSites.value);
        if (customRestrictionSites.length > 0) {
            payload.custom_restriction_sites = customRestrictionSites;
            state.customRestrictionSites = customRestrictionSites;
        }

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
                gc_bucket: gcBucket(primary.metrics.gc_percent ?? 0),
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
        elements.esmatlasFoldLink.href = `https://esmatlas.com/explore?tab=fold&sequence=${encoded}`;
    }
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
    elements.origCAI.textContent = 'N/A';
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

    // Sequence
    elements.optimizedSequence.textContent = formatSequence(primary.optimized_sequence);

    // Render GC Graph
    renderGCGraph(primary.optimized_sequence);
    renderCandidateComparison(res);
    renderCustomRestrictionResults(res);

    // PolyA color coding
    const polyaCount = primary.metrics.polya_signals;
    if (polyaCount === 0) {
        elements.polyaValue.textContent = '0 (Clean)';
        elements.polyaValue.className = 'text-xl font-black text-emerald-600 dark:text-emerald-400 leading-none';
    } else {
        elements.polyaValue.textContent = `${polyaCount} ⚠️`;
        elements.polyaValue.className = 'text-xl font-black text-amber-600 dark:text-amber-400 leading-none';
    }

    // Validation Badges
    elements.validationStatus.classList.remove('hidden');
    const v = primary.validation || { polya: 'PASS', moclo: 'UNCHECKED', gc: 'PASS' };
    updateValidationIcon('valPolyA', v.polya === 'PASS');
    updateValidationIcon('valMoClo', v.moclo === 'PASS' ? true : v.moclo === 'UNCHECKED' ? null : false);

    // Improved GC Status (Logic separation: N/A vs Out of Range)
    const valGC = document.getElementById('valGC');
    const gcMin = (primary.gc_window_min != null) ? Number(primary.gc_window_min) : 40.0;
    const gcMax = (primary.gc_window_max != null) ? Number(primary.gc_window_max) : 55.0;

    if (!optSeq || optSeq.length === 0 || calculatedGC === 0) {
        valGC.textContent = '⚠️';
        valGC.className = 'text-amber-500';
        valGC.nextElementSibling.textContent = 'GC Content: N/A';
        valGC.nextElementSibling.title = 'GC calculation failed or sequence not available';
    } else if (v.gc !== 'PASS' || (calculatedGC < gcMin || calculatedGC > gcMax)) {
        valGC.textContent = '⚠️';
        valGC.className = 'text-amber-500';
        valGC.nextElementSibling.textContent = `⚠️ Outside target range (${calculatedGC.toFixed(1)}%)`;
        valGC.nextElementSibling.title = `Target: ${gcMin.toFixed(0)}–${gcMax.toFixed(0)}% | Calculated: ${calculatedGC.toFixed(1)}%`;
    } else {
        valGC.textContent = '✅';
        valGC.className = 'text-emerald-400';
        valGC.nextElementSibling.textContent = `GC Content Check (${gcMin.toFixed(0)}–${gcMax.toFixed(0)}%)`;
    }

    // Structure prediction links
    updateStructureLinks();

    // JSON Details
    elements.jsonDetails.textContent = JSON.stringify(res, null, 2);
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
        aaPreserved: candidate.validator_status === 'pass' ? '✅ 100%' : '⚠️ Review'
    };
}

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
        const status = candidate.validator_status === 'pass' ? 'Pass' : 'Review';
        const aaStatus = candidate.validator_status === 'pass' ? '100%' : 'Review';
        return `
            <tr class="${isRecommended ? 'bg-emerald-50 dark:bg-emerald-900/20' : 'bg-white dark:bg-slate-900'}">
                <td class="px-3 py-2 font-bold text-slate-800 dark:text-slate-100">${candidate.label || candidate.id}${isRecommended ? ' ★' : ''}</td>
                <td class="px-3 py-2">${Number(candidate.cai || 0).toFixed(3)}</td>
                <td class="px-3 py-2">${Number(candidate.gc_percent || 0).toFixed(1)}</td>
                <td class="px-3 py-2">${Number(candidate.gc_window_min || 0).toFixed(1)}-${Number(candidate.gc_window_max || 0).toFixed(1)}</td>
                <td class="px-3 py-2">${aaStatus}</td>
                <td class="px-3 py-2">${candidate.internal_stop_count || 0}</td>
                <td class="px-3 py-2">${candidate.repeat_count || 0}</td>
                <td class="px-3 py-2 font-bold ${candidate.validator_status === 'pass' ? 'text-emerald-600' : 'text-amber-600'}">${status}</td>
            </tr>
        `;
    }).join('');
    elements.candidateComparisonContainer.classList.remove('hidden');
}

function renderCustomRestrictionResults(res) {
    if (!elements.customRestrictionResults || !elements.customRestrictionResultsBody) return;

    const custom = res.custom_restriction_sites;
    if (!custom) {
        elements.customRestrictionResults.classList.add('hidden');
        elements.customRestrictionResultsBody.innerHTML = '';
        return;
    }

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

function formatSequence(seq) {
    if (!seq) return '';
    return seq.match(/.{1,60}/g).join('\n');
}

function renderGCGraph(seq) {
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

    const ctx = elements.gcChart.getContext('2d');
    chartInstance = new Chart(ctx, {
        type: 'line',
        data: {
            labels: labels,
            datasets: [
                {
                    label: 'Target Max (65%)',
                    data: Array(n).fill(65),
                    borderColor: 'rgba(16, 185, 129, 0.35)',
                    borderWidth: 1,
                    borderDash: [4, 4],
                    pointRadius: 0,
                    fill: false,
                    tension: 0
                },
                {
                    label: 'Target Min (55%)',
                    data: Array(n).fill(55),
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
    };

    state.history.unshift(item);
    if (state.history.length > 10) state.history.pop();
    localStorage.setItem('factorforge_history', JSON.stringify(state.history));
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
        host_profile: item.host || 'nbenthamiana'
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
    state.sequence = '';
    state.customRestrictionSites = [];
    state.results = null;
    elements.previewContainer.classList.add('hidden');
    elements.inputTypeBadge.classList.add('hidden');
    updateInputStats('');
    elements.resultsContainer.classList.add('hidden');
    elements.constructIdDisplay.textContent = '';
    elements.constructIdRow.classList.add('hidden');
    if (elements.candidateComparisonContainer) elements.candidateComparisonContainer.classList.add('hidden');
    if (elements.customRestrictionResults) elements.customRestrictionResults.classList.add('hidden');
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
        const version = state.results.engine_versions?.product || '3.2.0';
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
