import { useEffect } from "react"
import { useMutation, useQuery } from "@apollo/react-hooks"
import {
  GET_ACTIVE_SPECTRA,
  GET_SELECTORS,
  NORMALIZE_CURRENT
} from "../client/queries"

// selectors is an array:
// [{id: 1, owner: String (slug from the ownersInfo dict), category: String}]
const useSelectors = ({ ownerInfo, spectraInfo }) => {
  const {
    data: { activeSpectra }
  } = useQuery(GET_ACTIVE_SPECTRA)
  const {
    data: { selectors }
  } = useQuery(GET_SELECTORS)

  const [normalize, { loading: normLoading }] = useMutation(NORMALIZE_CURRENT)
  useEffect(() => {
    if (!normLoading) normalize()
  }, [activeSpectra, normLoading, normalize])

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
