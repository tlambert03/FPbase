/**
 * Color palette definition for spectrum visualization.
 */
export interface Palette {
  /** Display name of the palette */
  name: string
  /** Array of hex color codes, or null for wavelength-based colors */
  hexlist: string[] | null
}

/**
 * Available color palettes for the spectra viewer.
 */
export type PaletteName = "wavelength" | "rainbow" | "tol_contrast" | "tol_vibrant" | "okabe_ito"

/**
 * Collection of predefined color palettes.
 */
export interface Palettes {
  wavelength: Palette
  // tableau20: Palette
  tableau20b: Palette
  rainbow: Palette
  tol_contrast: Palette
  tol_vibrant: Palette
  okabe_ito: Palette
}

const palettes: Palettes = {
  wavelength: {
    name: "By Wavelength",
    // this is the default palette, and uses colors based on wavelength from the database
    hexlist: null,
  },
  // tableau20: {
  //   name: "Tableau 20",
  //   hexlist: [
  //     "#1F77B4",
  //     "#AEC7E8",
  //     "#FF7F0E",
  //     "#FFBB78",
  //     "#2CA02C",
  //     "#98DF8A",
  //     "#D62728",
  //     "#FF9896",
  //     "#9467BD",
  //     "#C5B0D5",
  //     "#8C564B",
  //     "#C49C94",
  //     "#E377C2",
  //     "#F7B6D2",
  //     "#7F7F7F",
  //     "#C7C7C7",
  //     "#BCBD22",
  //     "#DBDB8D",
  //     "#17BECF",
  //     "#9EDAE5",
  //   ],
  // },
  tableau20b: {
    name: "Tableau 20B",
    hexlist: [
      "#393B79",
      "#5254A3",
      "#6B6ECF",
      "#9C9EDE",
      "#637939",
      "#8CA252",
      "#B5CF6B",
      "#CEDB9C",
      "#8C6D31",
      "#BD9E39",
      "#E7BA52",
      "#E7CB94",
      "#843C39",
      "#AD494A",
      "#D6616B",
      "#E7969C",
      "#7B4173",
      "#A55194",
      "#CE6DBD",
      "#DE9ED6",
    ],
  },
  rainbow: {
    name: "Rainbow (Paul Tol)",
    hexlist: [
      "#1965B0",
      "#DC050C",
      // "#F7F056", // COLOR 18
      "#F7CB45", // COLOR 19
      "#4EB265",
      "#7BAFDE",
      "#CAE0AB",
      "#882E72",
      "#EE8026",
      "#72190E",
      "#F4A736",
      "#5289c7",
      "#D1BBD7",
      "#90C987",
      "#AE76A3",
      "#F6C141",
      "#BA8DB4",
      "#994F88",
      "#A5170E",
    ],
  },
  // tol_muted: {
  //   name: 'Tol Muted',
  //   hexlist: [
  //     '#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933', '#CC6677',
  //     '#882255', '#AA4499', '#DDDDDD']
  // },
  // dark2: {
  //   name: 'Dark2 (matplotlib)',
  //   hexlist: [
  //     '#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d',
  //     '#666666']
  // },
  tol_contrast: {
    name: "High Contrast (Paul Tol; 4 color)",
    hexlist: ["#004488", "#DDAA33", "#BB5566", "#000000"],
  },
  tol_vibrant: {
    name: "Vibrant (Paul Tol)",
    hexlist: ["#0077BB", "#EE7733", "#33BBEE", "#CC3311", "#009988", "#EE3377", "#AAAAAA"],
  },
  // tol_light: {
  //   name: 'Tol Light',
  //   hexlist: [
  //     '#BBCC33', '#AAAA00', '#77AADD', '#EE8866', '#EEDD88', '#FFAABB', '#99DDFF',
  //     '#44BB99', '#DDDDDD']
  // },
  okabe_ito: {
    name: "Okabe Ito",
    hexlist: [
      "#E69F00",
      "#56B4E9",
      "#009E73",
      "#F0E442",
      "#0072B2",
      "#D55E00",
      "#CC79A7",
      "#000000",
    ],
  },
}

export default palettes
