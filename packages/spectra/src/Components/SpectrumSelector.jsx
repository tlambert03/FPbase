import Box from "@mui/material/Box"
import update from "immutability-helper"
import PropTypes from "prop-types"
import React, { useCallback, useEffect } from "react"
import { components } from "react-select"
import { useSpectrum } from "../hooks/useSpectraQueries"
import { useSpectraStore } from "../store/spectraStore"
import ProductLink from "./ProductLink"
import SortableWindowedSelect from "./SortableWindowedSelect"
import SubtypeSelector from "./SubtypeSelector"

import theme from "./theme"

const customStyles = {
  menu: (provided) => ({
    ...provided,
    zIndex: "10000",
  }),
  singleValue: (provided, _state) => ({
    ...provided,
    [theme.breakpoints.down("xs")]: {
      fontSize: "0.88rem",
    },
  }),
}

const SingleValue = ({ children, data, ...props }) => {
  const [extra, setExtra] = React.useState("")
  const spectrumId = data.category === "P" || data.category === "D" ? data.spectra[0]?.id : null
  const { data: spectrum } = useSpectrum(spectrumId)

  useEffect(() => {
    if (spectrum?.owner) {
      const { qy, extCoeff } = spectrum.owner
      if (qy || extCoeff) {
        let val = "("
        if (qy) {
          val += `QY: ${qy}`
          if (extCoeff) val += " / "
        }
        if (extCoeff) val += `EC: ${extCoeff.toLocaleString()}`
        val += ")"
        setExtra(val)
      } else {
        setExtra("")
      }
    }
  }, [spectrum])

  return (
    <components.SingleValue {...props}>
      {children}{" "}
      <span
        style={{
          fontSize: "0.76rem",
          color: "#a9a9a9",
          fontWeight: 600,
          marginLeft: 5,
          bottom: 1,
          position: "relative",
        }}
      >
        {extra}
      </span>
    </components.SingleValue>
  )
}
const SINGLE = { SingleValue }

const SpectrumSelector = React.memo(function SpectrumSelector({
  options,
  allOwners,
  showCategoryIcon,
  selector,
  ownerInfo,
}) {
  const [value, setValue] = React.useState({})

  useEffect(() => {
    setValue(selector.owner && ownerInfo[selector.owner])
  }, [ownerInfo, selector])

  const subtypes = value?.spectra || []

  // Get Zustand store actions and state
  const excludeSubtypes = useSpectraStore((state) => state.excludeSubtypes)
  const updateActiveSpectra = useSpectraStore((state) => state.updateActiveSpectra)

  // when the spectrum selector changes
  const handleOwnerChange = useCallback(
    (newValue) => {
      // if it's the same as the previous value do nothing
      if (value === newValue) return
      const newOwner = newValue?.value
      setValue(newOwner && ownerInfo[newOwner]) // FIXME: replace with optimistic UI?

      // Update active spectra - selectors will be derived automatically
      const spectraToAdd = newValue?.spectra
        .filter(({ subtype }) => !excludeSubtypes.includes(subtype))
        .map(({ id }) => id)
      const spectraToRemove = value?.spectra.map(({ id }) => id)

      updateActiveSpectra(spectraToAdd, spectraToRemove)
    },
    [excludeSubtypes, ownerInfo, value, updateActiveSpectra]
  )

  // disable options that are already claimed by other selectors

  const [myOptions, setMyOptions] = React.useState(options)

  useEffect(() => {
    const otherOwners = allOwners.filter((i) => i !== selector.owner)
    if (!otherOwners) return
    let newOptions = myOptions
    options.forEach((option, index) => {
      if (otherOwners.includes(option.value)) {
        if (!newOptions[index].isDisabled) {
          newOptions = update(newOptions, {
            [index]: {
              isDisabled: { $set: true },
              label: { $set: `${newOptions[index].label} (already selected)` },
            },
          })
        }
      } else if (newOptions[index].isDisabled) {
        newOptions = update(newOptions, {
          [index]: {
            isDisabled: { $set: false },
            label: { $set: option.label },
          },
        })
      }
    })
    if (newOptions !== myOptions) {
      setMyOptions(newOptions)
    }
  }, [allOwners, myOptions, options, selector.owner])

  return (
    <Box display="flex">
      <Box flexGrow={1}>
        <SortableWindowedSelect
          isClearable
          showIcon={showCategoryIcon}
          value={value}
          styles={customStyles}
          placeholder="Type to search..."
          onChange={handleOwnerChange}
          options={myOptions}
          components={SINGLE}
        />
      </Box>
      {subtypes.length > 0 && (
        <SubtypeSelector subtypes={subtypes} skip={value && !["P", "D"].includes(value.category)} />
      )}
      <ProductLink current={value} />
    </Box>
  )
})

SpectrumSelector.propTypes = {
  options: PropTypes.arrayOf(PropTypes.object).isRequired,
}

export default SpectrumSelector
