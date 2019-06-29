import { useEffect, useState } from "react"
import COLORS from "../../js/spectra/colors"
import { GET_SPECTRUM } from "../client/queries"
import { useApolloClient } from "react-apollo-hooks"

const rangexy = (start, end) =>
  Array.from({ length: end - start }, (v, k) => k + start)

// $cl1_wave
const customLaserSpectrum = _id => {
  let [id, wave] = _id.split("_")
  const data = [[+wave - 1, 0], [+wave, 1], [+wave + 1, 0]]
  const name = `${wave} laser`
  return Promise.resolve({
    data: {
      spectrum: {
        id: id,
        subtype: "L",
        owner: { name },
        category: "F",
        data,
        color: +wave in COLORS ? COLORS[+wave] : "#999999"
      }
    }
  })
}

// $cf1_type_center_width_trans
const customFilterSpectrum = _id => {
  let [id, subtype, center, width, trans] = _id.split("_")
  subtype = subtype.toUpperCase()
  trans = +trans / 100 || 0.9
  let data = []
  let name = `Custom `
  switch (subtype) {
    case "BP":
      const min = Math.round(+center - width / 2)
      const max = Math.round(+center + width / 2)
      data.push([min - 1, 0])
      rangexy(min, max + 1).forEach(x => data.push([x, +trans]))
      data.push([max + 1, 0])
      name += ` ${center}/${width} bp`
      break
    case "LP":
      rangexy(300, center).forEach(x => data.push([x, 0]))
      rangexy(+center + 1, 1000).forEach(x => data.push([x, +trans]))
      name += ` ${center}lp`
      break
    case "SP":
      rangexy(300, center).forEach(x => data.push([x, +trans]))
      rangexy(+center + 1, 1000).forEach(x => data.push([x, 0]))
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
        subtype,
        owner: { name },
        category: "F",
        data,
        color: +center in COLORS ? COLORS[+center] : "#999999"
      }
    }
  }
}

const useSpectralData = activeSpectra => {
  // $cf1_type_center_width
  const [data, setData] = useState([])
  const client = useApolloClient()
  useEffect(() => {
    async function fetchData() {
      let _data = await Promise.all(
        activeSpectra
          .filter(i => i)
          .map(id => {
            if (id.startsWith("$cf")) {
              return customFilterSpectrum(id)
            } else if (id.startsWith("$cl")) {
              return customLaserSpectrum(id)
            }
            return client.query({ query: GET_SPECTRUM, variables: { id } })
          })
      )
      setData(_data.map(({ data }) => data.spectrum).filter(i => i))
    }
    fetchData()
  }, [activeSpectra, client])

  return data
}

export default useSpectralData
