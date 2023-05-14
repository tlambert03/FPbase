import React, { useEffect, useCallback } from "react"
import PropTypes from "prop-types"
import Box from "@material-ui/core/Box"
import { useQuery, useMutation, useApolloClient } from "@apollo/react-hooks"
import { components } from "react-select"
import update from "immutability-helper"
import gql from "graphql-tag"
import { GET_OWNER_OPTIONS, GET_SPECTRUM } from "../client/queries"
import SubtypeSelector from "./SubtypeSelector"
import SortableWindowedSelect from "./SortableWindowedSelect"
import ProductLink from "./ProductLink"

import theme from "./theme"

const customStyles = {
  menu: provided => ({
    ...provided,
    zIndex: "10000",
  }),
  singleValue: (provided, state) => ({
    ...provided,
    [theme.breakpoints.down("xs")]: {
      fontSize: "0.88rem",
    },
  }),
}

const SingleValue = ({ children, data, ...props }) => {
  const client = useApolloClient()
  const [extra, setExtra] = React.useState("")
  useEffect(() => {
    async function fetchQeEc(id) {
      const {
        data: { spectrum },
      } = await client.query({
        query: GET_SPECTRUM,
        variables: { id: +id },
      })
      if (spectrum.owner) {
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
    }
    if (data.category === "P" || data.category === "D") {
      fetchQeEc(data.spectra[0].id)
    }
  }, [client, children, data.spectra, data.category])

  return (
    <components.SingleValue {...props}>
      {children}
      {" "}
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

const COMBO_MUTATE = gql`
  mutation comboMutate($add: [String], $remove: [String], $selector: Selector) {
    updateActiveSpectra(add: $add, remove: $remove) @client
    updateSelector(selector: $selector) @client
  }
`

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

  const subtypes = (value && value.spectra) || []
  // const [updateSpectra] = useMutation(UPDATE_ACTIVE_SPECTRA)

  const {
    data: { excludeSubtypes },
  } = useQuery(GET_OWNER_OPTIONS)

  // new Spectra get added in the handleOwnerChange handler in SpectrumSelector.jsx
  // const [updateSelector] = useMutation(UPDATE_SELECTOR)

  const [comboMutate] = useMutation(COMBO_MUTATE)
  // when the spectrum selector changes
  const handleOwnerChange = useCallback(
    newValue => {
      // if it's the same as the previous value do nothing
      if (value === newValue) return
      // setValue(newValue)
      // onChange(newValue && newValue.value)
      const newOwner = newValue && newValue.value
      setValue(newOwner && ownerInfo[newOwner]) // FIXME: replace with optimistic UI?

      comboMutate({
        variables: {
          add:
            newValue &&
            newValue.spectra
              .filter(({ subtype }) => !excludeSubtypes.includes(subtype))
              .map(({ id }) => id),
          remove: value && value.spectra.map(({ id }) => id),
          selector: {
            id: +selector.id,
            owner: newOwner,
            category: ownerInfo[newOwner] && ownerInfo[newOwner].category,
          },
        },
      })
    },
    [excludeSubtypes, ownerInfo, selector.id, value] // eslint-disable-line
  )

  // disable options that are already claimed by other selectors

  const [myOptions, setMyOptions] = React.useState(options)

  useEffect(() => {
    const otherOwners = allOwners.filter(i => i !== selector.owner)
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
        <SubtypeSelector
          subtypes={subtypes}
          skip={value && !["P", "D"].includes(value.category)}
        />
      )}
      <ProductLink current={value} />
    </Box>
  )
})

SpectrumSelector.propTypes = {
  options: PropTypes.arrayOf(PropTypes.object).isRequired,
}

export default SpectrumSelector
