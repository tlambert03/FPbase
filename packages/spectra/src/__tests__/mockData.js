import { List } from "immutable";

// Realistic spectrum data based on EGFP (excitation and emission)
export const mockSpectrumEGFP_EX = {
  id: 17,
  category: "P",
  subtype: "EX",
  color: "#00FF00",
  // Realistic EGFP excitation spectrum data: wavelength (nm), intensity (0-1)
  data: [
    [350, 0.01],
    [380, 0.05],
    [400, 0.15],
    [420, 0.35],
    [450, 0.75],
    [470, 0.95],
    [488, 1.0],
    [500, 0.85],
    [520, 0.25],
    [550, 0.05],
  ],
  owner: {
    id: "egfp-state",
    slug: "egfp",
    name: "EGFP",
    qy: 0.6,
    extCoeff: 55900,
    exMax: 488,
    emMax: 509,
  },
};

export const mockSpectrumEGFP_EM = {
  id: 18,
  category: "P",
  subtype: "EM",
  color: "#00FF00",
  // Realistic EGFP emission spectrum data
  data: [
    [470, 0.01],
    [490, 0.15],
    [500, 0.45],
    [509, 1.0],
    [520, 0.95],
    [540, 0.65],
    [560, 0.35],
    [580, 0.15],
    [600, 0.05],
    [620, 0.01],
  ],
  owner: {
    id: "egfp-state",
    slug: "egfp",
    name: "EGFP",
    qy: 0.6,
    extCoeff: 55900,
    exMax: 488,
    emMax: 509,
  },
};

export const mockChartOptions = {
  showY: false,
  showX: true,
  showGrid: false,
  areaFill: true,
  logScale: false,
  scaleEC: false,
  scaleQY: false,
  extremes: [null, null],
  shareTooltip: true,
  palette: "wavelength",
};

// Apollo Mock Provider setup
export const createMockApolloProvider = (spectraData = []) => {
  const mocks = {
    Query: {
      chartOptions: () => mockChartOptions,
      exNorm: () => [null],
      spectrum: (_, { id }) => {
        const spectrum = spectraData.find((s) => s.id === id);
        return spectrum || null;
      },
    },
  };
  return mocks;
};

// Convert regular array to Immutable List (simulates real app behavior)
export const toImmutableData = (data) => {
  return List(data.map((point) => List(point)));
};
