import gql from "graphql-tag"
import update from "immutability-helper"
import {
  GET_ACTIVE_SPECTRA,
  GET_EX_NORM,
  GET_SELECTORS,
  ADD_SELECTORS,
  batchSpectra,
  GET_ACTIVE_OVERLAPS,
} from "./queries"
import { trapz, spectraProduct } from "../util"
import PALETTES from "../palettes"

const PALETTE_KEYS = Object.keys(PALETTES)

export const defaults = {
  activeSpectra: [],
  activeOverlaps: [],
  chartOptions: {
    showY: false,
    showX: true,
    showGrid: false,
    logScale: false,
    scaleEC: false,
    scaleQY: false,
    shareTooltip: true,
    areaFill: true,
    palette: "wavelength",
    zoomType: "x",
    extremes: [undefined, undefined],
    __typename: "chartOptions",
  },
  exNorm: [null, null], // [normWave, normID]
  excludeSubtypes: ["2P"],
  selectors: [],
}

function toggleChartOption(cache, key) {
  const current = cache.readQuery({
    query: gql`
      {
        chartOptions {
          ${key}
        }
      }
    `,
  })
  const data = { ...current }
  data.chartOptions[key] = !current.chartOptions[key]
  cache.writeData({ data })
  return data
}

function _setPalette(palette, client) {
  const data = { chartOptions: { palette, __typename: "chartOptions" } }
  client.writeQuery({
    query: gql`
      {
        chartOptions @client {
          palette
        }
      }
    `,
    data,
  })
  return data
}

function activeSpectraToSelectors(
  activeSpectra,
  selectors,
  spectraInfo,
  ownerInfo
) {
  const currentOwners = selectors.map(({ owner }) => owner)
  const newOwners = activeSpectra
    .map(id => spectraInfo[id] && spectraInfo[id].owner)
    .filter(owner => owner && !currentOwners.includes(owner))
  const toAdd = [...new Set(newOwners)].map(owner => ({
    owner,
    category: ownerInfo[owner].category,
  }))
  Array.from(["D", "P", "L", "C", "F", null]).forEach(cat => {
    if (!selectors.find(item => item.category === cat && !item.owner)) {
      toAdd.push({
        owner: null,
        category: cat,
      })
    }
  })
  return toAdd
}

const isValidId = id => {
  if (!id) return false
  if (!Number.isNaN(parseFloat(id))) return true
  if (typeof id === "string") {
    if (id.startsWith("$cl") || id.startsWith("$cf")) {
      return true
    }
    return id.split("_").every(i => isValidId(i))
  }
  return false
}

const validSpectraIds = spectra => spectra.filter(id => isValidId(id))

export const resolvers = {
  Query: {
    overlap: async (_root, { ids }, { client }) => {
      const idString = ids.sort((a, b) => a - b).join("_")
      const { data } = await client.query({
        query: batchSpectra(ids),
      })
      const dataArray = ids.map(id => data[`spectrum_${id}`].data)
      const name = ids.map(id => data[`spectrum_${id}`].owner.name).join(" & ")
      const ownerID = ids
        .map(id => data[`spectrum_${id}`].owner.id)
        .sort((a, b) => a - b)
        .join("_")
      const product = spectraProduct(...dataArray)
      return {
        data: product,
        area: trapz(product),
        id: idString,
        category: "O",
        subtype: "O",
        color: "#000000",
        owner: { id: ownerID, name, __typename: "Owner" },
        __typename: "Spectrum",
      }
    },
  },
  Spectrum: {
    area: (spectrum, obj, cli) => {
      return trapz(spectrum.data)
    },
  },
  Mutation: {
    toggleYAxis: (_root, variables, { cache }) => {
      return toggleChartOption(cache, "showY")
    },
    toggleXAxis: (_root, variables, { cache }) => {
      return toggleChartOption(cache, "showX")
    },
    toggleGrid: (_root, variables, { cache }) => {
      return toggleChartOption(cache, "showGrid")
    },
    toggleLogScale: (_root, variables, { cache }) => {
      return toggleChartOption(cache, "logScale")
    },
    toggleScaleEC: (_root, variables, { cache }) => {
      return toggleChartOption(cache, "scaleEC")
    },
    toggleScaleQY: (_root, variables, { cache }) => {
      return toggleChartOption(cache, "scaleQY")
    },
    toggleShareTooltip: (_root, variables, { cache }) => {
      return toggleChartOption(cache, "shareTooltip")
    },
    toggleAreaFill: (_root, variables, { cache }) => {
      return toggleChartOption(cache, "areaFill")
    },
    setChartExtremes: (_root, { extremes }, { client }) => {
      const data = { chartOptions: { extremes, __typename: "chartOptions" } }
      client.writeQuery({
        query: gql`
          {
            chartOptions @client {
              extremes
            }
          }
        `,
        data,
      })
      return data
    },
    setPalette: (_root, { palette }, { client }) => {
      return _setPalette(palette, client)
    },
    cyclePalette: (_root, variables, { cache, client }) => {
      const {
        chartOptions: { palette },
      } = cache.readQuery({
        query: gql`
          {
            chartOptions {
              palette
            }
          }
        `,
      })
      const curIndex = PALETTE_KEYS.indexOf(palette)
      const newpalette = PALETTE_KEYS[(curIndex + 1) % PALETTE_KEYS.length]
      return _setPalette(newpalette, client)
    },
    setExcludeSubtypes: (_, { excludeSubtypes }, { cache }) => {
      cache.writeData({ data: { excludeSubtypes } })
    },
    setExNorm: async (_, { data }, { client }) => {
      await client.writeQuery({ query: GET_EX_NORM, data: { exNorm: data } })
      return data
    },
    setActiveSpectra: async (_, { activeSpectra }, { cache, client }) => {
      const filtered = [...new Set(activeSpectra)]
        // .filter(id => Boolean(spectrumFrag(cache, id)))
        .map(i => String(i))
      const data = {
        activeSpectra: validSpectraIds(filtered),
      }
      await client.writeQuery({ query: GET_ACTIVE_SPECTRA, data })
      return data
    },
    updateActiveSpectra: async (_, { add, remove }, { cache, client }) => {
      let { activeSpectra } = cache.readQuery({ query: GET_ACTIVE_SPECTRA })
      activeSpectra = activeSpectra.filter(id => {
        if (id.startsWith("$cf") && remove) {
          const _id = id.split("_")[0]
          return !(remove.findIndex(item => item.startsWith(_id)) > -1)
        }
        if (id.startsWith("$cl") && remove) {
          const _id = id.split("_")[0]
          return !(remove.findIndex(item => item.startsWith(_id)) > -1)
        }
        return !(remove || []).includes(id)
      })
      const toAdd = (add || []).filter(id => id).map(id => String(id))
      const data = {
        activeSpectra: validSpectraIds([
          ...new Set([...activeSpectra, ...toAdd]),
        ]),
      }
      await client.writeQuery({ query: GET_ACTIVE_SPECTRA, data })
      // client.mutate({
      //   mutation: NORMALIZE_CURRENT
      // })
      return data
    },
    updateActiveOverlaps: async (_, { add, remove }, { cache, client }) => {
      let { activeOverlaps } = cache.readQuery({ query: GET_ACTIVE_OVERLAPS })
      activeOverlaps = activeOverlaps.filter(id => !(remove || []).includes(id))
      const toAdd = (add || []).filter(id => id)
      const data = {
        activeOverlaps: validSpectraIds([
          ...new Set([...activeOverlaps, ...toAdd]),
        ]),
      }
      await client.writeQuery({ query: GET_ACTIVE_OVERLAPS, data })
      return data
    },
    normalizeCurrent: (_, args, { cache, client }) => {
      const { activeSpectra, selectors: currentSelectors } = cache.readQuery({
        query: gql`
          {
            activeSpectra @client
            selectors @client
          }
        `,
      })
      const selectors = activeSpectraToSelectors(
        activeSpectra,
        currentSelectors,
        window.spectraInfo,
        window.ownerInfo
      )
      if (selectors.length > 0) {
        client.mutate({
          mutation: ADD_SELECTORS,
          variables: { selectors },
        })
      }
      return [currentSelectors, selectors]
    },
    addSelectors: (_, { selectors }, { cache, client }) => {
      const { selectors: currentSelectors } = cache.readQuery({
        query: GET_SELECTORS,
      })
      let selectorIDs = currentSelectors.reduce(
        (acc, next) => Math.max(acc, next.id),
        0
      )
      const newSelectors = selectors.map(sel => {
        sel.id = ++selectorIDs
        return sel
      })
      const data = {
        selectors: update(currentSelectors, { $push: newSelectors }),
      }
      client.writeQuery({ query: GET_SELECTORS, data })
    },
    updateSelector: (_, { selector }, { cache, client }) => {
      const { selectors } = cache.readQuery({ query: GET_SELECTORS })
      const index = selectors.findIndex(item => item.id === selector.id)
      let data
      if (selector.owner) {
        data = {
          selectors: update(selectors, { [index]: { $set: selector } }),
        }
      } else {
        data = {
          selectors: update(selectors, { $splice: [[index, 1]] }),
        }
      }
      client.writeQuery({ query: GET_SELECTORS, data })
    },
    removeSelector: (_, { id }, { cache, client }) => {
      const { selectors } = cache.readQuery({ query: GET_SELECTORS })
      const index = selectors.findIndex(selector => selector.id === id)
      const data = {
        selectors: update(selectors, { $splice: [[index, 1]] }),
      }
      client.writeQuery({ query: GET_SELECTORS, data })
    },
    clearForm: async (_, args, { cache }) => {
      const { leave, appendSpectra } = args || {}
      const { activeSpectra } = cache.readQuery({ query: GET_ACTIVE_SPECTRA })
      let keepSpectra = []
      if ((leave || []).length > 0) {
        keepSpectra = activeSpectra.filter(
          id =>
            window.spectraInfo[id] &&
            leave.includes(window.spectraInfo[id].category)
        )
      }
      cache.writeData({
        data: {
          exNorm: [null, null],
          selectors: [],
          activeOverlaps: [],
          activeSpectra: [
            ...new Set([...keepSpectra, ...(appendSpectra || [])]),
          ],
        },
      })
    },
  },
}

export { validSpectraIds }

// function spectrumFrag(cache, id) {
//   return cache.readFragment({
//     id: "SpectrumInfo:" + String(id),
//     fragment: gql`
//       fragment Spectrum on SpectrumInfo {
//         owner {
//           slug
//         }
//       }
//     `
//   })
// }
