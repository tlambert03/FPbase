import gql from "graphql-tag"
import { GET_ACTIVE_SPECTRA, GET_EX_NORM } from "./queries"

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
    extremes: [null, null],
    __typename: "chartOptions"
  },
  exNorm: [null, null], // [normWave, normID]
  excludeSubtypes: ["2P"]
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

export const resolvers = {
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
    setActiveSpectra: async (_, { activeSpectra }, { cache, client }) => {
      const filtered = [...new Set(activeSpectra)]
        //.filter(id => Boolean(spectrumFrag(cache, id)))
        .map(i => String(i))
      const data = {
        activeSpectra: filtered
      }
      await client.writeQuery({ query: GET_ACTIVE_SPECTRA, data })
      return data
    },
    setExcludeSubtypes: (_, { excludeSubtypes }, { cache }) => {
      cache.writeData({ data: { excludeSubtypes } })
    },
    updateActiveSpectra: async (_, { add, remove }, { cache, client }) => {
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
        activeSpectra: [...new Set([...activeSpectra, ...toAdd])]
      }
      await client.writeQuery({ query: GET_ACTIVE_SPECTRA, data })
      return data
    },
    setExNorm: async (_, { data }, { client }) => {
      await client.writeQuery({ query: GET_EX_NORM, data: { exNorm: data } })
      return data
    }
  }
}
