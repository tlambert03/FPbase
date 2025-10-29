import { useApolloClient, useQuery } from "@apollo/client"
import gql from "graphql-tag"
import update from "immutability-helper"
import { useEffect, useState } from "react"
import { GET_SPECTRUM } from "../client/queries"
import COLORS from "../colors"

const rangexy = (start, end) => Array.from({ length: end - start }, (_v, k) => k + start)

// $cl1_wave
const customLaserSpectrum = (_id) => {
  let [id, wave] = _id.split("_")

  wave = +wave
  const data = [
    [wave - 1, 0],
    [wave, 1],
    [wave + 1, 0],
  ]
  const name = `${wave} laser`
  return Promise.resolve({
    data: {
      spectrum: {
        id: id,
        customId: _id,
        subtype: "PD",
        owner: { name, id: _id, slug: _id },
        category: "L",
        data,
        color: wave in COLORS ? COLORS[wave] : "#999999",
      },
    },
  })
}

// $cf1_type_center_width_trans
const customFilterSpectrum = (_id) => {
  let [id, subtype, center, width, trans] = _id.split("_")

  subtype = subtype.toUpperCase()
  trans = +trans / 100 || 0.9
  const data = []
  let name = `Custom `
  switch (subtype) {
    case "BP": {
      const min = Math.round(+center - width / 2)
      const max = Math.round(+center + width / 2)
      data.push([min - 1, 0])
      rangexy(min, max + 1).forEach((x) => {
        data.push([x, +trans])
      })
      data.push([max + 1, 0])
      name += ` ${center}/${width} bp`
      break
    }
    case "LP":
      rangexy(300, center).forEach((x) => {
        data.push([x, 0])
      })
      rangexy(+center + 1, 1000).forEach((x) => {
        data.push([x, +trans])
      })
      name += ` ${center}lp`
      break
    case "SP":
      rangexy(300, center).forEach((x) => {
        data.push([x, +trans])
      })
      rangexy(+center + 1, 1000).forEach((x) => {
        data.push([x, 0])
      })
      name += ` ${center}sp`
      break
    default:
      break
  }

  return {
    data: {
      spectrum: {
        // setting this to "id" is faster, but causes an error when you mash the sliders
        // setting it to "_id" is safer, but incurs a full update with each change
        id: id,
        customId: _id,
        subtype,
        owner: { name, id: _id, slug: _id },
        category: "F",
        data,
        color: (+center) in COLORS ? COLORS[+center] : "#999999",
      },
    },
  }
}

const useSpectralData = (provideSpectra, provideOverlaps) => {
  // $cf1_type_center_width
  const [currentData, setCurrentData] = useState([])
  const client = useApolloClient()
  const { data } = useQuery(
    gql`
      {
        activeSpectra @client
        activeOverlaps @client
      }
    `,
    { skip: provideSpectra }
  )
  let activeSpectra
  let activeOverlaps
  if (provideSpectra) {
    activeSpectra = provideSpectra
    activeOverlaps = provideOverlaps
  } else {
    activeSpectra = data?.activeSpectra || []
    activeOverlaps = data?.activeOverlaps || []
  }
  useEffect(() => {
    function idToData(id) {
      // cast id to integer
      if (id.startsWith("$cf")) {
        return customFilterSpectrum(id)
      }
      if (id.startsWith("$cl")) {
        return customLaserSpectrum(id)
      }
      return client.query({ query: GET_SPECTRUM, variables: { id: +id } })
    }

    async function updateData() {
      // Use functional state update to get current data without having it as a dependency
      // This prevents the infinite loop risk
      setCurrentData((prevData) => {
        // Compute changes based on previous state
        const deadSpectra = prevData.reduceRight((acc, item, idx) => {
          if (
            !activeSpectra.includes(item.customId || item.id) &&
            !activeOverlaps.includes(item.id)
          )
            acc.push([idx, 1])
          return acc
        }, [])

        const currentIDs = prevData.map((item) => item.customId || item.id)
        const newSpectraIds = activeSpectra.filter((id) => id && !currentIDs.includes(id))
        const newOverlapIds = activeOverlaps.filter((id) => id && !currentIDs.includes(id))

        // If no changes, return the same reference (no re-render)
        if (!deadSpectra.length && !newSpectraIds.length && !newOverlapIds.length) {
          return prevData
        }

        // Schedule async fetch for new data
        if (newSpectraIds.length || newOverlapIds.length) {
          ;(async () => {
            let newData = await Promise.all(newSpectraIds.map((id) => idToData(id)))
            newData = newData.map((item) => item.data.spectrum).filter((i) => i)

            const newOverlapData = newOverlapIds
              .map((id) => window.OverlapCache[id])
              .filter((i) => i)

            // Add fetched data
            if (newData.length || newOverlapData.length) {
              setCurrentData((current) =>
                update(current, {
                  $push: [...newData, ...newOverlapData],
                })
              )
            }
          })()
        }

        // Immediately remove dead spectra
        if (deadSpectra.length) {
          return update(prevData, { $splice: deadSpectra })
        }

        return prevData
      })
    }

    updateData()
  }, [activeOverlaps, activeSpectra, client]) // Removed currentData to prevent infinite loop

  return currentData
}

export default useSpectralData
