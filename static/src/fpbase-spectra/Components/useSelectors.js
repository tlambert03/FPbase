import { useState, useRef, useEffect } from "react"
import { useMutation, useQuery } from "react-apollo-hooks"
import {
  UPDATE_ACTIVE_SPECTRA,
  GET_ACTIVE_SPECTRA,
  SET_ACTIVE_SPECTRA
} from "../client/queries"
import update from "immutability-helper"

const useSelectors = ({ owners, spectraInfo, initial = [] }) => {
  const {
    data: { activeSpectra }
  } = useQuery(GET_ACTIVE_SPECTRA)
  const [selectors, setSelectors] = useState(initial)
  const selectorId = useRef(0)

  // this makes sure that all active spectra are reflected in the selectors array
  useEffect(() => {
    const currentOwners = selectors.map(({ owner }) => owner)
    const newOwners = activeSpectra
      .map(id => spectraInfo[id] && spectraInfo[id].owner)
      .filter(owner => owner && !currentOwners.includes(owner))
    let newSelectors = [...new Set(newOwners)].map(owner => ({
      id: selectorId.current++,
      owner,
      category: owners[owner].category
    }))

    newSelectors = update(selectors, { $push: newSelectors })

    setSelectors(newSelectors)
  }, [activeSpectra, owners, selectors, spectraInfo])

  useEffect(() => {
    let newSelectors = selectors
    Array.from(["D", "P", "L", "C", "F", null]).forEach(cat => {
      if (!newSelectors.find(item => item.category === cat && !item.owner)) {
        newSelectors = update(newSelectors, {
          $push: [
            {
              id: selectorId.current++,
              owner: null,
              category: cat
            }
          ]
        })
      }
    })
    setSelectors(newSelectors)
  }, [selectors])

  const addRow = (category = null) => {
    let newSelectors = update(selectors, {
      $push: [
        {
          id: selectorId.current++,
          owner: null,
          category
        }
      ]
    })
    setSelectors(newSelectors)
  }

  const updateSpectra = useMutation(UPDATE_ACTIVE_SPECTRA)
  const removeRow = selector => {
    if (owners[selector.owner] && owners[selector.owner].spectra) {
      updateSpectra({
        variables: {
          remove: owners[selector.owner].spectra.map(({ id }) => id)
        }
      })
    }
    setSelectors(selectors.filter(({ id }) => id !== selector.id))
  }

  const changeOwner = (id, category = null) => {
    return function(newOwner) {
      const index = selectors.findIndex(selector => selector.id === id)
      let newSelectors
      if (newOwner) {
        newSelectors = update(selectors, {
          [index]: {
            owner: { $set: newOwner },
            category: {
              $set: newOwner ? owners[newOwner].category : category
            }
          }
        })
        setSelectors(newSelectors)
      } else {
        removeRow({ id, category, owner: newOwner })
      }
    }
  }

  const setSpectra = useMutation(SET_ACTIVE_SPECTRA)
  const clearForm = (leave = [], appendSpectra = []) => {
    const preserve = selectors.filter(({ category }) =>
      leave.includes(category)
    )
    setSelectors([
      ...preserve,
      {
        id: selectorId.current++,
        owner: null,
        category: null
      }
    ])
    const keepSpectra = activeSpectra.filter(
      id => spectraInfo[id] && leave.includes(spectraInfo[id].category)
    )
    setSpectra({
      variables: {
        activeSpectra: [...new Set([...keepSpectra, ...appendSpectra])]
      }
    })
  }

  return { selectors, addRow, clearForm, changeOwner, removeRow }
}

export default useSelectors
