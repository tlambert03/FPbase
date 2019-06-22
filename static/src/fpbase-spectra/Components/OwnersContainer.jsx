import React, { useState, useRef, useEffect } from "react"
import PropTypes from "prop-types"
import { useMutation, useQuery } from "react-apollo-hooks"
import {
  UPDATE_ACTIVE_SPECTRA,
  GET_ACTIVE_SPECTRA,
  SET_ACTIVE_SPECTRA
} from "../client/queries"
import update from "immutability-helper"
import Tabs from "@material-ui/core/Tabs"
import Tab from "@material-ui/core/Tab"
import { makeStyles } from "@material-ui/core/styles"
import MyAppBar from "./MyAppBar"

import SpectrumSelectorGroup from "./SpectrumSelectorGroup"
import { Typography } from "@material-ui/core"
import useWindowWidth from "./useWindowWidth"

const useStyles = makeStyles(theme => ({
  tabHeader: {
    //marginBottom: 12,
    //marginLeft: 60,
    //paddingLeft: 0
    whiteSpace: "nowrap"
  },
  tabLabel: {
    marginTop: 0,
    paddingTop: 0,
    minHeight: 40,
    // lineHeight: 0,
    [theme.breakpoints.down("xs")]: {
      fontSize: "0.7rem",
      paddingLeft: 5,
      paddingRight: 5
    },
    [theme.breakpoints.up("sm")]: {
      fontSize: ".73rem",
      paddingLeft: "4%",
      paddingRight: "4%"
    },
    [theme.breakpoints.up("md")]: {
      fontSize: ".76rem",
      paddingLeft: 20,
      paddingRight: 20
    },
    [theme.breakpoints.up("lg")]: {
      fontSize: ".9rem",
      paddingLeft: 60,
      paddingRight: 60
    }
  },
  categoryHeader: {
    textTransform: "uppercase",
    fontSize: "small",
    color: "#3F51B5",
    marginTop: ".2rem",
    marginBottom: "0.4rem"
  }
}))

function selectorSorter(a, b) {
  const ORDER = ["P", "D", "F", "L", "C", "", null]
  if (ORDER.indexOf(a.category) === ORDER.indexOf(b.category)) {
    if (a.owner && !b.owner) return -1
    if (!a.owner && b.owner) return 1
    return 0
  }
  if (ORDER.indexOf(a.category) > ORDER.indexOf(b.category)) return 1
  return -1
}

const OwnersContainer = ({ owners, spectraInfo }) => {
  const {
    data: { activeSpectra }
  } = useQuery(GET_ACTIVE_SPECTRA)
  const [selectors, setSelectors] = useState([])
  const selectorId = useRef(0)
  const classes = useStyles()

  const [tab, setTab] = useState(0)

  const handleTabChange = (event, newValue) => {
    if (newValue !== tab) {
      setTab(newValue)
    }
  }

  // keyboard shortcut for tab switcher
  useEffect(() => {
    const handleKeyDown = event => {
      // don't do anything if we're on an input
      if (document.activeElement.tagName.toUpperCase() === "INPUT") {
        return
      }
      switch (event.code) {
        case "Digit1":
        case "Digit2":
        case "Digit3":
        case "Digit4":
        case "Digit5":
          setTab(event.key - 1)
          break
        case "ArrowRight":
          setTab(prevTab => (prevTab === 4 ? 0 : prevTab + 1))
          break
        case "ArrowLeft":
          setTab(prevTab => (prevTab === 0 ? 4 : prevTab - 1))
          break
        default:
          break
      }
    }

    document.addEventListener("keydown", handleKeyDown)
    return () => {
      document.removeEventListener("keydown", handleKeyDown)
    }
  }, [])

  // this makes sure that all active spectra are reflected in the ownerSlugs array
  useEffect(() => {
    const currentOwners = selectors.map(({ owner }) => owner)
    let newOwners = activeSpectra
      .map(id => spectraInfo[id] && spectraInfo[id].owner)
      .filter(owner => owner && !currentOwners.includes(owner))
    newOwners = [...new Set(newOwners)].map(owner => ({
      id: selectorId.current++,
      owner,
      category: owners[owner].category
    }))
    setSelectors(update(selectors, { $push: newOwners }))
  }, [activeSpectra]) // eslint-disable-line

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

  const isPopulated = cat =>
    selectors.filter(({ owner, category }) => category === cat && owner)
      .length > 0

  const smartLabel = (label, cats) => {
    let populated
    if (cats === null) {
      populated = Boolean(selectors.filter(i => i.owner).length)
    } else {
      populated = cats.some(c => isPopulated(c))
    }
    return (
      <span
        style={{
          fontWeight: populated ? "bold" : "normal",
          whiteSpace: "nowrap"
        }}
      >
        {label}
        {populated ? " ✶" : ""}
      </span>
    )
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

  selectors.sort(selectorSorter)
  const allOptions = Object.values(owners)
  const spectrumCategoryGroup = (category, hint) => {
    return (
      <SpectrumSelectorGroup
        selectors={selectors}
        options={allOptions}
        addRow={addRow}
        showCategoryIcon={!Boolean(category)}
        changeOwner={changeOwner}
        removeRow={removeRow}
        owners={owners}
        category={category}
        hint={hint}
      />
    )
  }

  const width = useWindowWidth()

  return (
    <div style={{ paddingBottom: 100 }}>
      <Tabs
        value={tab}
        onChange={handleTabChange}
        indicatorColor="primary"
        textColor="primary"
        centered={width >= 500 ? true : false}
        variant={width >= 500 ? "standard" : "scrollable"}
        scrollButtons="on"
        className={classes.tabHeader}
      >
        <Tab           tabIndex={-1} className={classes.tabLabel} label={smartLabel("All", null)} />
        <Tab
          className={classes.tabLabel}
          tabIndex={-1}
          label={smartLabel("Fluorophores", ["D", "P"])}
        />
        <Tab
          className={classes.tabLabel}
          tabIndex={-1}
          label={smartLabel("Filters", ["F"])}
        />
        <Tab
          tabIndex={-1}
          className={classes.tabLabel}
          label={smartLabel("Lights", ["L"])}
        />
        <Tab
          tabIndex={-1}
          className={classes.tabLabel}
          label={smartLabel("Detectors", ["C"])}
        />
      </Tabs>

      <TabContainer index={tab}>
        <div>{spectrumCategoryGroup()}</div>
        <div>
          <Typography variant="h6" className={classes.categoryHeader}>
            Fluorescent Proteins
          </Typography>
          {spectrumCategoryGroup("P", "protein")}
          <Typography variant="h6" className={classes.categoryHeader}>
            Dyes
          </Typography>
          {spectrumCategoryGroup("D", "dye")}
        </div>
        <div>{spectrumCategoryGroup("F", "filter")}</div>
        <div>{spectrumCategoryGroup("L", "light")}</div>
        <div>{spectrumCategoryGroup("C", "camera")}</div>
      </TabContainer>

      <MyAppBar spectraOptions={allOptions} clearForm={clearForm} />
    </div>
  )
}

OwnersContainer.propTypes = {
  category: PropTypes.string
}

OwnersContainer.defaultProps = {
  category: ""
}

// here to make sure we always render each tab to maintain the formstate
const TabContainer = ({ index, children }) =>
  children && React.cloneElement(children[index])

export default OwnersContainer
