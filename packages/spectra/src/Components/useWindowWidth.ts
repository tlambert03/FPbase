import { useEffect, useState } from "react"
import { debounce } from "../util"

/**
 * Custom hook that tracks the current window width with debounced resize events.
 *
 * @param _delay - Unused parameter (kept for backwards compatibility)
 * @returns The current window width in pixels
 */
const useWindowWidth = (_delay?: number): number => {
  const [width, setWidth] = useState<number>(window.innerWidth)

  useEffect(() => {
    const handleResize = debounce(() => {
      setWidth(window.innerWidth)
    }, 100)

    window.addEventListener("resize", handleResize)
    return () => {
      window.removeEventListener("resize", handleResize)
    }
  }, []) // Empty dependency array - only run on mount/unmount

  return width
}

export default useWindowWidth
