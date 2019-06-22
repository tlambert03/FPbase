import React, { useEffect } from "react"
import PropTypes from "prop-types"
import Box from "@material-ui/core/Box"
import { useQuery, useMutation, useApolloClient } from "react-apollo-hooks"
import {
  UPDATE_ACTIVE_SPECTRA,
  GET_OWNER_OPTIONS,
  GET_SPECTRUM
} from "../client/queries"
import SubtypeSelector from "./SubtypeSelector"
import SortablePaginatedSelect from "./SortablePaginatedSelect"
import ProductLink from "./ProductLink"
import { components } from "react-select"

import theme from "./theme"
const customStyles = {
  control: provided => ({
    ...provided,
    [theme.breakpoints.down("sm")]: {
      minHeight: 34,
      height: 34
    }
  }),
  menu: provided => ({
    ...provided,
    zIndex: "10000"
  }),
  placeholder: provided => ({
    ...provided,
    [theme.breakpoints.down("sm")]: {
      top: "42%",
      fontSize: "0.9rem"
    }
  }),
  indicatorsContainer: provided => ({
    ...provided,
    [theme.breakpoints.down("sm")]: {
      position: "relative",
      top: "-1px",
      height: 34
    }
  }),
  singleValue: (provided, state) => ({
    ...provided,
    [theme.breakpoints.down("sm")]: {
      fontSize: "0.9rem"
    }
  })
}

const SingleValue = ({ children, ...props }) => {
  const client = useApolloClient()
  const [extra, setExtra] = React.useState("")
  useEffect(() => {
    async function fetchQeEc(id) {
      const { data } = await client.query({
        query: GET_SPECTRUM,
        variables: { id }
      })
      if (data.spectrum.owner) {
        const { qy, extCoeff } = data.spectrum.owner
        if (qy || extCoeff) {
          let val = "("
          if (qy) {
            val += `QY: ${qy}`
            if (extCoeff) val += " / "
          }
          if (extCoeff) val += `EC: ${extCoeff.toLocaleString()}`
          val += ")"
          setExtra(val)
        }
      }
    }
    if (props.data.category === "P" || props.data.category === "D") {
      fetchQeEc(props.data.spectra[0].id)
    }
  }, [client, children, props.data.spectra, props.data.category])

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
          position: "relative"
        }}
      >
        {extra}
      </span>
    </components.SingleValue>
  )
}

const SpectrumSelector = ({
  options,
  value,
  otherOwners,
  onChange,
  showCategoryIcon
}) => {
  //const [value, setValue] = useState(current)
  const subtypes = (value && value.spectra) || []
  const updateSpectra = useMutation(UPDATE_ACTIVE_SPECTRA)

  const {
    data: { excludeSubtypes }
  } = useQuery(GET_OWNER_OPTIONS)

  // when the spectrum selector changes
  const handleOwnerChange = newValue => {
    // if it's the same as the previous value do nothing
    if (value === newValue) return
    // setValue(newValue)
    onChange(newValue && newValue.value)
    updateSpectra({
      variables: {
        add:
          newValue &&
          newValue.spectra
            .filter(({ subtype }) => !excludeSubtypes.includes(subtype))
            .map(({ id }) => id),
        remove: value && value.spectra.map(({ id }) => id)
      }
    })
  }

  // disable options that are already claimed by other selectors
  const myOptions = options.map(opt =>
    (otherOwners || []).includes(opt.value)
      ? {
          ...opt,
          label: opt.label + " (already selected)",
          isDisabled: true
        }
      : opt
  )

  return (
    <Box display="flex">
      <Box flexGrow={1}>
        <SortablePaginatedSelect
          isClearable
          showIcon={showCategoryIcon}
          value={value}
          styles={customStyles}
          placeholder="Type to search..."
          onChange={handleOwnerChange}
          options={myOptions}
          components={{ SingleValue }}
        />
      </Box>
      <SubtypeSelector
        subtypes={subtypes}
        skip={value && !["P", "D"].includes(value.category)}
      />
      <ProductLink current={value} />
    </Box>
  )
}

SpectrumSelector.propTypes = {
  options: PropTypes.arrayOf(PropTypes.object).isRequired,
  current: PropTypes.objectOf(PropTypes.any),
  otherOwners: PropTypes.arrayOf(PropTypes.string),
  showIcon: PropTypes.bool
}

SpectrumSelector.defaultProps = {
  otherOwners: [],
  current: null,
  showIcon: true
}

export default SpectrumSelector
