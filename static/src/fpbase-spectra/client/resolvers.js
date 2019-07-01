import gql from "graphql-tag"
import {
  GET_ACTIVE_SPECTRA,
  GET_EX_NORM,
  GET_SELECTORS,
  ADD_SELECTORS
} from "./queries"
import update from "immutability-helper"

export const defaults = {
  activeSpectra: [],
  chartOptions: {
    showY: false,
    showX: true,
    showGrid: false,
    logScale: false,
    scaleEC: false,
    scaleQY: false,
    shareTooltip: true,
    areaFill: true,
    extremes: [undefined, undefined],
    __typename: "chartOptions"
  },
  exNorm: [null, null], // [normWave, normID]
  excludeSubtypes: ["2P"],
  selectors: []
}

function toggleChartOption(cache, key) {
  const current = cache.readQuery({
    query: gql`
      {
        chartOptions {
          ${key}
        }
      }
    `
  })
  const data = { ...current }
  data.chartOptions[key] = !current.chartOptions[key]
  cache.writeData({ data })
  return data
}

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

function trapz(arr) {
  // approximate area under curve as series of trapezoids
  // arr = [[wave, data], ...]
  let sum = 0
  for (let i = 1; i < arr.length; i++) {
    sum += 0.5 * (arr[i][1] + arr[i - 1][1]) * (arr[i][0] - arr[i - 1][0])
  }
  return sum
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
  let toAdd = [...new Set(newOwners)].map(owner => ({
    owner,
    category: ownerInfo[owner].category
  }))
  Array.from(["D", "P", "L", "C", "F", null]).forEach(cat => {
    if (!selectors.find(item => item.category === cat && !item.owner)) {
      toAdd.push({
        owner: null,
        category: cat
      })
    }
  })
  return toAdd
}

export const resolvers = {
  Spectrum: {
    area: (spectrum, obj, cli) => {
      return trapz(spectrum.data)
    }
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
        data
      })
      return data
    },
    setExcludeSubtypes: (_, { excludeSubtypes }, { cache }) => {
      cache.writeData({ data: { excludeSubtypes } })
    },
    setExNorm: async (_, { data }, { client }) => {
      await client.writeQuery({ query: GET_EX_NORM, data: { exNorm: data } })
      return data
    },
    setActiveSpectra: async (_, { activeSpectra }, { cache, client }) => {
      console.timeLog("timer", "setActiveSpectra")
      const filtered = [...new Set(activeSpectra)]
        //.filter(id => Boolean(spectrumFrag(cache, id)))
        .map(i => String(i))
      const data = {
        activeSpectra: validSpectraIds(filtered)
      }
      await client.writeQuery({ query: GET_ACTIVE_SPECTRA, data })
      return data
    },
    updateActiveSpectra: async (_, { add, remove }, { cache, client }) => {
      console.timeLog("timer", "updateActiveSpectra")
      let { activeSpectra } = cache.readQuery({ query: GET_ACTIVE_SPECTRA })
      activeSpectra = activeSpectra.filter(id => {
        if (id.startsWith("$cf") && remove) {
          const _id = id.split("_")[0]
          return !(remove.findIndex(item => item.startsWith(_id)) > -1)
        } else if (id.startsWith("$cl") && remove) {
          const _id = id.split("_")[0]
          return !(remove.findIndex(item => item.startsWith(_id)) > -1)
        }
        return !(remove || []).includes(id)
      })
      const toAdd = (add || []).filter(id => id).map(id => String(id))
      const data = {
        activeSpectra: validSpectraIds([
          ...new Set([...activeSpectra, ...toAdd])
        ])
      }
      await client.writeQuery({ query: GET_ACTIVE_SPECTRA, data })
      // client.mutate({
      //   mutation: NORMALIZE_CURRENT
      // })
      return data
    },
    normalizeCurrent: (_, args, { cache, client }) => {
      console.timeLog("timer", "normalize Called")
      let { activeSpectra, selectors: currentSelectors } = cache.readQuery({
        query: gql`
          {
            activeSpectra @client
            selectors @client
          }
        `
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
          variables: { selectors }
        })
      }
      return [currentSelectors, selectors]
    },
    addSelectors: (_, { selectors }, { cache, client }) => {
      let { selectors: currentSelectors } = cache.readQuery({
        query: GET_SELECTORS
      })
      let selectorIDs = currentSelectors.reduce(
        (acc, next) => Math.max(acc, next.id),
        0
      )
      selectors = selectors.map(sel => {
        sel.id = ++selectorIDs
        return sel
      })
      const data = { selectors: update(currentSelectors, { $push: selectors }) }
      console.timeLog("timer", "addSelectors", data.selectors)
      client.writeQuery({ query: GET_SELECTORS, data })
    },
    updateSelector: (_, { selector }, { cache, client }) => {
      console.timeLog("timer", "updateSelectors")
      let { selectors } = cache.readQuery({ query: GET_SELECTORS })
      const index = selectors.findIndex(item => item.id === selector.id)
      let data
      if (selector.owner) {
        data = {
          selectors: update(selectors, { [index]: { $set: selector } })
        }
      } else {
        data = {
          selectors: update(selectors, { $splice: [[index, 1]] })
        }
      }
      client.writeQuery({ query: GET_SELECTORS, data })
    },
    removeSelector: (_, { id }, { cache, client }) => {
      console.timeLog("timer", "removeSelector")
      let { selectors } = cache.readQuery({ query: GET_SELECTORS })
      const index = selectors.findIndex(selector => selector.id === id)
      const data = {
        selectors: update(selectors, { $splice: [[index, 1]] })
      }
      client.writeQuery({ query: GET_SELECTORS, data })
    }
  }
}

const validSpectraIds = spectra =>
  spectra.filter(
    id => id && (!isNaN(id) || id.startsWith("$cl") || id.startsWith("$cf"))
  )

export { validSpectraIds }
