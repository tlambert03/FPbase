import { useMutation, useQuery } from "@apollo/client"
import { useEffect, useState } from "react"
import { GET_ACTIVE_SPECTRA, SET_ACTIVE_SPECTRA } from "../client/queries"

const CurrentSpectraInput = () => {
  const { loading, data } = useQuery(GET_ACTIVE_SPECTRA)
  const activeSpectra = data?.activeSpectra || []

  const [updateSpectra] = useMutation(SET_ACTIVE_SPECTRA)
  const [value, setValue] = useState("")
  useEffect(() => {
    setValue(activeSpectra.join(", "))
  }, [activeSpectra])
  const handleKeyPress = (e) => {
    if (e.key !== "Enter") {
      return
    }
    e.preventDefault()
    const newVal = e.target.value
      .split(",")
      .map((i) => i.trim())
      .filter((i) => i)
    updateSpectra({
      variables: {
        activeSpectra: newVal,
      },
      update: (
        _cache,
        {
          data: {
            setActiveSpectra: { activeSpectra },
          },
        }
      ) => {
        setValue(activeSpectra.join(", "))
      },
    })
  }
  const handleChange = (e) => {
    setValue(e.target.value)
  }
  if (loading) {
    return <></>
  }
  return (
    <input
      type="text"
      value={value}
      style={{ width: "80%" }}
      onKeyPress={handleKeyPress}
      onChange={handleChange}
    />
  )
}

export default CurrentSpectraInput
