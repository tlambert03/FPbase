import { useEffect, useState } from "react"
import { useSpectraStore } from "../store/spectraStore"

const CurrentSpectraInput = () => {
  const activeSpectra = useSpectraStore((state) => state.activeSpectra)
  const updateActiveSpectra = useSpectraStore((state) => state.updateActiveSpectra)

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
    updateActiveSpectra(newVal)
    setValue(newVal.join(", "))
  }
  const handleChange = (e) => {
    setValue(e.target.value)
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
