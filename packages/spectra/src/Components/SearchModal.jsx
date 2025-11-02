import Checkbox from "@mui/material/Checkbox"
import FormControlLabel from "@mui/material/FormControlLabel"
import Modal from "@mui/material/Modal"
import Typography from "@mui/material/Typography"
import { makeStyles } from "@mui/styles"
import PropTypes from "prop-types"
import React, { useEffect, useState } from "react"
import { components } from "react-select"
import { fetchGraphQL } from "../api/client"
import { GET_OPTICAL_CONFIG } from "../api/queries"
import { useSpectraStore } from "../store/spectraStore"
import { useCachedFetch } from "../useCachedQuery"
import MuiReactSelect from "./MuiReactSelect"

/**
 * Filter spectra by category
 * @param {string[]} spectraIds - Array of spectrum IDs
 * @param {Object} spectraInfo - Mapping of spectrum ID to spectrum metadata
 * @param {string[]} keepCategories - Categories to keep (e.g., ["P", "D"])
 * @returns {string[]} Filtered spectrum IDs
 */
function filterSpectraByCategory(spectraIds, spectraInfo, keepCategories) {
  if (!spectraInfo) return []
  return spectraIds.filter((id) => {
    const spectrum = spectraInfo[id]
    return spectrum && keepCategories.includes(spectrum.category)
  })
}

const OptionWithBlurb = (props) => {
  const myProps = { ...props }

  myProps.children = (
    <>
      {myProps.children}
      <span style={{ fontSize: "0.69rem", color: "#666" }}>
        <br />
        {myProps.data.comments}
      </span>
    </>
  )
  return <components.Option {...myProps} />
}

function getModalStyle() {
  const top = 40
  const left = 42

  return {
    top: `${top}%`,
    left: `${left}%`,
    transform: `translate(-${top}%, -${left}%)`,
  }
}

const useStyles = makeStyles((theme) => ({
  paper: {
    position: "absolute",
    maxWidth: "620px",
    width: "620px",
    backgroundColor: theme.palette.background.paper,
    boxShadow: theme.shadows[5],
    padding: "5px 30px 15px",
    outline: "none",
    [theme.breakpoints.down("sm")]: {
      width: "80%",
    },
    [theme.breakpoints.down("xs")]: {
      width: "84%",
      padding: "5px 22px 15px",
      fontSize: "small",
    },
  },
  title: {
    color: "black",
    marginTop: 29,
    marginBotton: 5,
    [theme.breakpoints.down("xs")]: {
      fontSize: "1.2rem",
    },
  },
}))

const filterOption = ({ label }, query) => {
  if (!query) return true
  const words = query.toLowerCase().split(" ")
  const opt = label.toLowerCase()
  return words.reduce((acc, cur) => acc && opt.includes(cur), true)
}

const SearchModal = React.memo(function SearchModal({ options, open, setOpen }) {
  const [modalStyle] = useState(getModalStyle)
  const classes = useStyles()

  const [ocOptions, setOcOptions] = useState([])
  const stash = useCachedFetch("/api/proteins/ocinfo/", "_FPbaseOCStash", 10 * 60)
  // const stash = useCachedQuery(OPTICAL_CONFIG_LIST, "_FPbaseOCStash", 5 * 60)
  useEffect(() => {
    if (stash) {
      const newOpts = stash.opticalConfigs.map(({ id, name, microscope, comments }) => ({
        label: `${name} (${microscope.name})`,
        value: id,
        comments: comments,
      }))
      setOcOptions(newOpts)
    }
  }, [stash])

  useEffect(() => {
    const handleKeyDown = (event) => {
      // don't do anything if we're on an input
      if (document.activeElement && document.activeElement.tagName.toUpperCase() === "INPUT") {
        return
      }
      switch (event.code) {
        // case "KeyL":
        case "Space":
          event.preventDefault()
          setOpen(true)
          break
        case "Escape":
          event.preventDefault()
          setOpen(false)
          break
        default:
          break
      }
    }

    document.addEventListener("keydown", handleKeyDown)
    return () => {
      document.removeEventListener("keydown", handleKeyDown)
    }
  }, [setOpen])

  const excludeSubtypes = useSpectraStore((state) => state.excludeSubtypes)
  const activeSpectra = useSpectraStore((state) => state.activeSpectra)
  const setActiveSpectra = useSpectraStore((state) => state.setActiveSpectra)
  const updateActiveSpectra = useSpectraStore((state) => state.updateActiveSpectra)

  const handleChange = async (event) => {
    const spectra = event?.spectra
      .filter(({ subtype }) => !excludeSubtypes.includes(subtype))
      .map(({ id }) => id)

    if (spectra && event) {
      // Just update activeSpectra - selectors will be derived automatically
      updateActiveSpectra(spectra)
    }
    setOpen(false)
  }
  const [preserveFluors, setPreserveFluors] = useState(true)
  const handleOCChange = async ({ value }) => {
    const { opticalConfig } = await fetchGraphQL(GET_OPTICAL_CONFIG, { id: value })
    const newSpectra = []
    opticalConfig.filters &&
      newSpectra.push(...opticalConfig.filters.map(({ spectrum }) => spectrum.id))
    opticalConfig.light && newSpectra.push(opticalConfig.light.spectrum.id)
    opticalConfig.camera && newSpectra.push(opticalConfig.camera.spectrum.id)

    // Keep fluorophores if preserveFluors is true, otherwise clear all
    const keptSpectra = preserveFluors
      ? filterSpectraByCategory(activeSpectra, window.spectraInfo, ["P", "D"])
      : []

    setActiveSpectra([...keptSpectra, ...newSpectra])
    setOpen(false)
  }

  return (
    <Modal open={open} onClose={() => setOpen(false)}>
      <div style={modalStyle} className={classes.paper}>
        <Typography variant="h6" className={classes.title}>
          Quick Entry
        </Typography>

        <MuiReactSelect
          autoFocus
          showIcon
          options={options}
          onChange={handleChange}
          className={classes.select}
          placeholder="Search for any spectrum"
        />

        <Typography variant="h6" className={classes.title}>
          Optical Config Lookup
        </Typography>
        <MuiReactSelect
          paginate={false}
          isLoading={ocOptions.length === 0}
          options={ocOptions}
          onChange={handleOCChange}
          className={classes.select}
          filterOption={filterOption}
          components={{ Option: OptionWithBlurb }}
          placeholder="Search optical configs"
        />
        <p
          style={{
            fontSize: "smaller",
            color: "#666",
            fontStyle: "italic",
            marginTop: 8,
            marginBottom: 1,
          }}
        >
          This is a list of optical configurations used in FPbase user microscopes. Selecting one
          will populate the chart with items from that config.
        </p>
        <FormControlLabel
          control={
            <Checkbox
              onChange={() => setPreserveFluors(!preserveFluors)}
              checked={preserveFluors}
              value="preserveFluors"
            />
          }
          label={<span style={{ fontSize: "small" }}>Preserve fluorophores</span>}
        />
        <p
          style={{
            marginTop: 8,
            fontSize: "small",
            textAlign: "center",
          }}
        >
          Spectra missing?{" "}
          <a href="/spectra/submit/" target="_blank" rel="noopener">
            Submit a spectrum
          </a>
        </p>
      </div>
    </Modal>
  )
})

SearchModal.propTypes = {
  options: PropTypes.arrayOf(PropTypes.object).isRequired,
}

export default SearchModal
