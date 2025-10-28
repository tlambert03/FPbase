/**
 * Test fixtures for spectrum data
 * These mirror real data structures from the FPbase GraphQL API
 */

export const EGFP_STATE_OWNER = {
  __typename: 'State',
  id: 'egfp-state',
  slug: 'egfp',
  name: 'EGFP',
  // Fields from FluorophoreInterface fragment
  qy: 0.6,
  extCoeff: 55900,
  twopPeakgm: null,
  exMax: 488,
  emMax: 509,
}

export const MCHERRY_STATE_OWNER = {
  __typename: 'State',
  id: 'mcherry-state',
  slug: 'mcherry',
  name: 'mCherry',
  // Fields from FluorophoreInterface fragment
  qy: 0.22,
  extCoeff: 72000,
  twopPeakgm: null,
  exMax: 587,
  emMax: 610,
}

export const ALEXA488_DYE_OWNER = {
  __typename: 'Dye',
  id: 'alexa-488',
  slug: 'alexa-488',
  name: 'Alexa Fluor 488',
  // Fields from FluorophoreInterface fragment
  qy: 0.92,
  extCoeff: 71000,
  twopPeakgm: null,
  exMax: 495,
  emMax: 519,
}

export const FILTER_OWNER = {
  __typename: 'Filter',
  id: 'filter-1',
  slug: 'bp-525-50',
  name: 'BP 525/50',
  // Filters don't implement FluorophoreInterface
  // so they don't have qy, extCoeff, etc.
}

export const EGFP_EXCITATION_SPECTRUM = {
  __typename: 'Spectrum',
  id: 18,
  category: 'P',
  subtype: 'EX',
  color: '#00FF00',
  data: [
    [350, 0.01],
    [400, 0.15],
    [450, 0.45],
    [488, 1.0],
    [520, 0.35],
    [550, 0.05],
  ],
  owner: EGFP_STATE_OWNER,
}

export const EGFP_EMISSION_SPECTRUM = {
  __typename: 'Spectrum',
  id: 17,
  category: 'P',
  subtype: 'EM',
  color: '#00FF00',
  data: [
    [450, 0.01],
    [480, 0.05],
    [509, 1.0],
    [550, 0.45],
    [600, 0.05],
  ],
  owner: EGFP_STATE_OWNER,
}

export const MCHERRY_EXCITATION_SPECTRUM = {
  __typename: 'Spectrum',
  id: 80,
  category: 'P',
  subtype: 'EX',
  color: '#FF0066',
  data: [
    [450, 0.05],
    [500, 0.15],
    [550, 0.55],
    [587, 1.0],
    [620, 0.25],
  ],
  owner: MCHERRY_STATE_OWNER,
}

export const MCHERRY_EMISSION_SPECTRUM = {
  __typename: 'Spectrum',
  id: 79,
  category: 'P',
  subtype: 'EM',
  color: '#FF0066',
  data: [
    [550, 0.01],
    [580, 0.08],
    [610, 1.0],
    [650, 0.45],
    [700, 0.05],
  ],
  owner: MCHERRY_STATE_OWNER,
}

export const FILTER_SPECTRUM = {
  __typename: 'Spectrum',
  id: 50,
  category: 'F',
  subtype: 'BP',
  color: '#0000FF',
  data: [
    [500, 0.0],
    [510, 0.9],
    [530, 0.9],
    [540, 0.0],
  ],
  owner: FILTER_OWNER,
}

/**
 * Helper to get spectrum by ID (mimics API behavior)
 */
export function getSpectrumById(id) {
  const spectra = {
    17: EGFP_EMISSION_SPECTRUM,
    18: EGFP_EXCITATION_SPECTRUM,
    79: MCHERRY_EMISSION_SPECTRUM,
    80: MCHERRY_EXCITATION_SPECTRUM,
    50: FILTER_SPECTRUM,
  }
  return spectra[id] || null
}
