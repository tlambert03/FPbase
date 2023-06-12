import { useEffect, useState } from "react"
import { useApolloClient, useQuery } from "@apollo/react-hooks"
import update from "immutability-helper"
import gql from "graphql-tag"
import { GET_SPECTRUM } from "../client/queries"
import COLORS from "../colors"

const rangexy = (start, end) =>
  Array.from({ length: end - start }, (v, k) => k + start)

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
        owner: { name, id: _id },
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
      rangexy(min, max + 1).forEach((x) => data.push([x, +trans]))
      data.push([max + 1, 0])
      name += ` ${center}/${width} bp`
      break
    }
    case "LP":
      rangexy(300, center).forEach((x) => data.push([x, 0]))
      rangexy(+center + 1, 1000).forEach((x) => data.push([x, +trans]))
      name += ` ${center}lp`
      break
    case "SP":
      rangexy(300, center).forEach((x) => data.push([x, +trans]))
      rangexy(+center + 1, 1000).forEach((x) => data.push([x, 0]))
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
        owner: { name, id: _id },
        category: "F",
        data,
        color: +center in COLORS ? COLORS[+center] : "#999999",
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
    activeSpectra = data.activeSpectra
    activeOverlaps = data.activeOverlaps
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
      // find dead Spectra
      const deadSpectra = currentData.reduceRight((acc, item, idx) => {
        if (
          !activeSpectra.includes(item.customId || item.id) &&
          !activeOverlaps.includes(item.id)
        )
          acc.push([idx, 1])
        return acc
      }, [])

      // find new activeSpectra that aren't in current Data
      const currentIDs = currentData.map((item) => item.customId || item.id)
      const newSpectra = activeSpectra.filter(
        (id) => id && !currentIDs.includes(id)
      )
      let newData = await Promise.all(newSpectra.map((id) => idToData(id)))
      newData = newData.map((item) => item.data.spectrum).filter((i) => i)
      // find new overlaps that aren't in current Data
      const newOverlaps = activeOverlaps.filter(
        (id) => id && !currentIDs.includes(id)
      )
      const newOverlapData = newOverlaps
        .map((id) => window.OverlapCache[id] || id)
        .filter((i) => i)

      if (deadSpectra.length || newData.length || newOverlapData.length) {
        setCurrentData(
          update(currentData, {
            $splice: deadSpectra,
            $push: [...newData, ...newOverlapData],
          })
        )
      }
    }

    updateData()
  }, [activeOverlaps, activeSpectra, client, currentData])

  return currentData
}

export default useSpectralData
