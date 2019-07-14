import { useEffect } from "react"
import { useMutation } from "@apollo/react-hooks"
import gql from "graphql-tag"

const toggleMut = param => gql`
mutation Toggle${param} {
  toggle${param} @client
}
`

const useKeyboardShortcuts = () => {
  const [toggleY] = useMutation(toggleMut("YAxis"))
  const [toggleX] = useMutation(toggleMut("XAxis"))
  const [toggleGrid] = useMutation(toggleMut("Grid"))
  const [toggleScaleEC] = useMutation(toggleMut("ScaleEC"))
  const [toggleScaleQY] = useMutation(toggleMut("ScaleQY"))
  const [toggleShareTooltip] = useMutation(toggleMut("ShareTooltip"))
  const [toggleAreaFill] = useMutation(toggleMut("AreaFill"))

  useEffect(() => {
    const handleKeyDown = event => {
      // don't do anything if we're on an input
      if (
        document.activeElement &&
        document.activeElement.tagName.toUpperCase() === "INPUT"
      ) {
        return
      }
      switch (event.code) {
        case "KeyX":
          event.preventDefault()
          toggleX()
          break
        case "KeyY":
          event.preventDefault()
          toggleY()
          break
        case "KeyG":
          event.preventDefault()
          toggleGrid()
          break
        case "KeyE":
          event.preventDefault()
          toggleScaleEC()
          break
        case "KeyQ":
          event.preventDefault()
          toggleScaleQY()
          break
        case "KeyA":
        case "KeyF":
          event.preventDefault()
          toggleAreaFill()
          break
        case "KeyT":
          event.preventDefault()
          toggleShareTooltip()
          break
        default:
          break
      }
    }

    document.addEventListener("keydown", handleKeyDown)
    return () => {
      document.removeEventListener("keydown", handleKeyDown)
    }
  }, [
    toggleAreaFill,
    toggleGrid,
    toggleScaleEC,
    toggleScaleQY,
    toggleShareTooltip,
    toggleX,
    toggleY,
  ])
}

export default useKeyboardShortcuts
