import { useEffect, useState } from 'react'
import { debounce } from '../util'

const useWindowWidth = (_delay) => {
  const [width, setWidth] = useState(window.innerWidth)

  useEffect(() => {
    const handleResize = debounce(() => {
      setWidth(window.innerWidth)
    }, 100)

    window.addEventListener('resize', handleResize)
    return () => {
      window.removeEventListener('resize', handleResize)
    }
  }, []) // Empty dependency array - only run on mount/unmount

  return width
}

export default useWindowWidth
