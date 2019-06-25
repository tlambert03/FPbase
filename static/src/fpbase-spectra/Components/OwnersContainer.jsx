import React, { useState, useEffect } from "react"
import PropTypes from "prop-types"
import Tabs from "@material-ui/core/Tabs"
import Tab from "@material-ui/core/Tab"
import { makeStyles } from "@material-ui/core/styles"
import { RingLoader } from "react-spinners"
import { css } from "@emotion/core"
import SpectrumSelectorGroup from "./SpectrumSelectorGroup"
import { Typography } from "@material-ui/core"
import useWindowWidth from "./useWindowWidth"
import CustomFilterGroup from "./CustomFilterGroup"
import CustomLaserGroup from "./CustomLaserGroup"

const override = css`
  display: block;
  margin: 50px auto;
  border-color: red;
`

const useStyles = makeStyles(theme => ({
  tabHeader: {
    //marginBottom: 12,
    //marginLeft: 60,
    //paddingLeft: 0
    whiteSpace: "nowrap",
    minHeight: 52
  },
  tabLabel: {
    marginTop: 0,
    paddingTop: 3,
    paddingBottom: 3,
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
    marginTop: ".35rem",
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

const OwnersContainer = ({
  owners,
  selectors,
  changeOwner,
  removeRow,
  activeSpectra
}) => {
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

  const isPopulated = cat =>
    selectors.filter(({ owner, category }) => category === cat && owner)
      .length > 0

  const smartLabel = (label, cats) => {
    let populated
    if (cats === null) {
      populated = Boolean(selectors.filter(i => i.owner).length)
    } else {
      populated = cats.some(c => isPopulated(c))
      if (cats.find(i => i === "F")) {
        populated = populated || activeSpectra.some(s => s.startsWith("$cf"))
      }
      if (cats.find(i => i === "L")) {
        populated = populated || activeSpectra.some(s => s.startsWith("$cl"))
      }
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

  selectors.sort(selectorSorter)
  const allOptions = Object.values(owners)
  const spectrumCategoryGroup = (category, hint) => {
    return (
      <SpectrumSelectorGroup
        selectors={selectors}
        options={allOptions}
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
        <Tab
          tabIndex={-1}
          className={classes.tabLabel}
          label={smartLabel("All", null)}
        />
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
          label={smartLabel("Light Sources", ["L"])}
        />
        <Tab
          tabIndex={-1}
          className={classes.tabLabel}
          label={smartLabel("Detectors", ["C"])}
        />
      </Tabs>

      {Object.keys(owners).length === 0 ? (
        <div className="sweet-loading">
          <RingLoader
            css={override}
            sizeUnit={"px"}
            size={100}
            color={"#ccc"}
            loading={true}
          />
        </div>
      ) : (
        <TabContainer index={tab}>
          <div>
            {spectrumCategoryGroup()}
          </div>
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
          <div>
            {spectrumCategoryGroup("F", "filter")}
            <CustomFilterGroup activeSpectra={activeSpectra} />
          </div>
          <div>
            {spectrumCategoryGroup("L", "light")}
            <CustomLaserGroup activeSpectra={activeSpectra} />
          </div>
          <div>{spectrumCategoryGroup("C", "camera")}</div>
        </TabContainer>
      )}
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
