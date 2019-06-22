import React, { useState, useEffect } from "react"
import PropTypes from "prop-types"
import Modal from "@material-ui/core/Modal"
import Typography from "@material-ui/core/Typography"
import { makeStyles } from "@material-ui/core/styles"
import MuiReactSelect from "./MuiReactSelect"
import { useMutation, useQuery, useApolloClient } from "react-apollo-hooks"
import { components } from "react-select"
import {
  UPDATE_ACTIVE_SPECTRA,
  GET_OWNER_OPTIONS,
  OPTICAL_CONFIG_LIST,
  GET_OPTICAL_CONFIG
} from "../client/queries"
import { useCachedQuery } from "../useCachedQuery"
import Checkbox from "@material-ui/core/Checkbox"
import FormControlLabel from "@material-ui/core/FormControlLabel"

function getModalStyle() {
  const top = 40
  const left = 42

  return {
    top: `${top}%`,
    left: `${left}%`,
    transform: `translate(-${top}%, -${left}%)`
  }
}

const useStyles = makeStyles(theme => ({
  paper: {
    position: "absolute",
    maxWidth: "620px",
    width: "620px",
    backgroundColor: theme.palette.background.paper,
    boxShadow: theme.shadows[5],
    padding: 30,
    paddingTop: 5,
    outline: "none",
    [theme.breakpoints.down("sm")]: {
      width: "80%"
    },
    [theme.breakpoints.down("xs")]: {
      width: "84%",
      padding: "5px 22px 30px",
      fontSize: "small"
    }
  },
  title: {
    color: "black",
    marginTop: 25,
    marginBotton: 5,
    [theme.breakpoints.down("xs")]: {
      fontSize: "1.2rem"
    }
  }
}))

const filterOption = ({ label }, query) => {
  if (!query) return true
  const words = query.toLowerCase().split(" ")
  const opt = label.toLowerCase()
  return words.reduce((acc, cur) => acc && opt.includes(cur), true)
}

const SearchModal = ({ options, open, setOpen, clearForm }) => {
  const [modalStyle] = useState(getModalStyle)
  const classes = useStyles()

  const [ocOptions, setOcOptions] = useState([])
  const stash = useCachedQuery(OPTICAL_CONFIG_LIST, "_FPbaseOCStash", 5 * 60)
  useEffect(() => {
    if (stash) {
      const newOpts = stash.opticalConfigs.map(
        ({ id, name, microscope, comments }) => ({
          label: `${name} (${microscope.name})`,
          value: id,
          comments: comments
        })
      )
      setOcOptions(newOpts)
    }
  }, [stash])

  useEffect(() => {
    const handleKeyDown = event => {
      // don't do anything if we're on an input
      if (document.activeElement.tagName.toUpperCase() === "INPUT") {
        return
      }
      switch (event.code) {
        case "KeyL":
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

  const updateSpectra = useMutation(UPDATE_ACTIVE_SPECTRA)

  const {
    data: { excludeSubtypes }
  } = useQuery(GET_OWNER_OPTIONS)

  const handleChange = event => {
    const spectra =
      event &&
      event.spectra
        .filter(({ subtype }) => !excludeSubtypes.includes(subtype))
        .map(({ id }) => id)

    if (spectra) {
      updateSpectra({ variables: { add: spectra } })
    }
    setOpen(false)
  }

  const client = useApolloClient()
  const [preserveFluors, setPreserveFluors] = useState(true)
  const handleOCChange = async ({ value }) => {
    const {
      data: { opticalConfig }
    } = await client.query({
      query: GET_OPTICAL_CONFIG,
      variables: { id: value }
    })
    const spectra = []
    opticalConfig.filters &&
      spectra.push(...opticalConfig.filters.map(({ spectrum }) => spectrum.id))
    opticalConfig.light && spectra.push(opticalConfig.light.spectrum.id)
    opticalConfig.camera && spectra.push(opticalConfig.camera.spectrum.id)

    clearForm(preserveFluors ? ["P", "D"] : [], spectra)
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
        <p style={{ fontSize: "smaller", color: "#666", fontStyle: "italic" }}>
          This is a list of optical configurations used in FPbase user
          microcopes. Selecting one will load all items from that config into
          the spectra viewer.
        </p>
        <FormControlLabel
          control={
            <Checkbox
              onChange={() => setPreserveFluors(!preserveFluors)}
              checked={preserveFluors}
              value={"preserveFluors"}
            />
          }
          label="Preserve fluorophores when changing config"
        />
      </div>
    </Modal>
  )
}

const OptionWithBlurb = props => {
  const myProps = { ...props }

  myProps.children = (
    <>
      {myProps.children}
      <span style={{ fontSize: "0.69rem", color: "#666" }}>
        <br />
        {props.data.comments}
      </span>
    </>
  )
  return <components.Option {...myProps} />
}

SearchModal.propTypes = {
  options: PropTypes.arrayOf(PropTypes.object).isRequired
}

export default SearchModal
