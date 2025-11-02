import { useEffect } from "react"
import { useSpectraStore } from "../store/spectraStore"

// selectors is an array:
// [{id: 1, owner: String (slug from the ownersInfo dict), category: String}]
const useSelectors = () => {
  const activeSpectra = useSpectraStore((state) => state.activeSpectra)
  const selectors = useSpectraStore((state) => state.selectors)
  const normalizeCurrent = useSpectraStore((state) => state.normalizeCurrent)

  useEffect(() => {
    normalizeCurrent()
  }, [normalizeCurrent])

  return { activeSpectra, selectors }
}

export default useSelectors

// const [setSpectra] = useMutation(SET_ACTIVE_SPECTRA)
// const [mutateExNormWave] = useMutation(SET_EX_NORM)
// const clearForm = useCallback(
//   (leave = [], appendSpectra = []) => {
//     const preserve = selectors.filter(({ category }) =>
//       leave.includes(category)
//     )
//     setSelectors([
//       ...preserve,
//       {
//         id: selectorId.current++,
//         owner: null,
//         category: null
//       }
//     ])
//     const keepSpectra = activeSpectra.filter(
//       id => spectraInfo[id] && leave.includes(spectraInfo[id].category)
//     )
//     mutateExNormWave({ variables: { data: [null, null] } })
//     setSpectra({
//       variables: {
//         activeSpectra: [...new Set([...keepSpectra, ...appendSpectra])]
//       }
//     })
//   },
//   [activeSpectra, mutateExNormWave, selectors, setSpectra, spectraInfo]
// )
