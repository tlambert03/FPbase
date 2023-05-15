import React, { useState, useEffect, useMemo } from "react"
import Tabs from "@material-ui/core/Tabs"
import Tab from "@material-ui/core/Tab"
import { makeStyles } from "@material-ui/core/styles"
import { RingLoader } from "react-spinners"
import { css } from "@emotion/core"
import { Typography } from "@material-ui/core"
import { useQuery, useApolloClient } from "@apollo/react-hooks"
import gql from "graphql-tag"
import SpectrumSelectorGroup from "./SpectrumSelectorGroup"
import CustomFilterGroup from "./CustomFilterGroup"
import CustomLaserGroup from "./CustomLaserGroup"
// import useSelectors from "./useSelectors"
import { NORMALIZE_CURRENT } from "../client/queries"
import { isTouchDevice } from "../util"
import EfficiencyTable from "./EfficiencyTable"
import { categoryIcon } from "./FaIcon"

const ISTOUCH = isTouchDevice()

const override = css`
  display: block;
  margin: 50px auto;
  border-color: red;
`

const useStyles = makeStyles(theme => ({
  tabHeader: {
    marginBottom: 12,
    // marginLeft: 60,
    // paddingLeft: 0
    whiteSpace: "nowrap",
    minHeight: 40,
  },
  tabLabel: {
    marginTop: 0,
    paddingTop: 3,
    paddingBottom: 3,
    minHeight: 40,
    minWidth: 72,
    // lineHeight: 0,
    [theme.breakpoints.down("xs")]: {
      fontSize: "0.7rem",
      paddingLeft: 5,
      paddingRight: 5,
    },
    [theme.breakpoints.up("sm")]: {
      fontSize: ".73rem",
      paddingLeft: "4%",
      paddingRight: "4%",
    },
    [theme.breakpoints.up("md")]: {
      fontSize: ".76rem",
      paddingLeft: 20,
      paddingRight: 20,
      minWidth: 160,
    },
    [theme.breakpoints.up("lg")]: {
      fontSize: ".9rem",
      paddingLeft: 60,
      paddingRight: 60,
    },
  },
  bigShow: {
    [theme.breakpoints.down("sm")]: {
      display: "none",
    },
    [theme.breakpoints.up("md")]: {
      display: "block",
    },
  },
  bigHide: {
    [theme.breakpoints.down("sm")]: {
      display: "block",
      marginRight: 13,
      marginLeft: 13,
    },
    [theme.breakpoints.down("xs")]: {
      display: "block",
    },
    [theme.breakpoints.up("md")]: {
      display: "none",
    },
  },
  categoryHeader: {
    textTransform: "uppercase",
    fontSize: "small",
    color: "#3F51B5",
    marginTop: ".35rem",
    marginBottom: "0.4rem",
  },
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

const OwnersContainer = React.memo(function OwnersContainer({
  ownerInfo,
  spectraInfo,
}) {
  const classes = useStyles()
  const [tab, setTab] = useState(0)

  const {
    data: { activeSpectra, selectors },
  } = useQuery(gql`
    {
      selectors @client
      activeSpectra @client
    }
  `)

  const client = useApolloClient()
  useEffect(() => {
    client.mutate({ mutation: NORMALIZE_CURRENT })
  }, [activeSpectra, client]) //eslint-disable-line

  const handleTabChange = (event, newValue) => {
    if (newValue !== tab) {
      setTab(newValue)
    }
  }

  // keyboard shortcut for tab switcher
  useEffect(() => {
    const handleKeyDown = event => {
      // don't do anything if we're on an input
      if (
        document.activeElement &&
        document.activeElement.tagName.toUpperCase() === "INPUT"
      ) {
        return
      }
      switch (event.code) {
        case "Digit1":
        case "Digit2":
        case "Digit3":
        case "Digit4":
        case "Digit5":
        case "Digit6":
          setTab(event.key - 1)
          break
        case "ArrowRight":
          setTab(prevTab => (prevTab === 5 ? 0 : prevTab + 1))
          break
        case "ArrowLeft":
          setTab(prevTab => (prevTab === 0 ? 5 : prevTab - 1))
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

  const isPopulated = cat => {
    let populated =
      selectors.filter(({ owner, category }) => category === cat && owner)
        .length > 0
    if (cat === "F") {
      populated = populated || activeSpectra.some(s => s.startsWith("$cf"))
    }
    if (cat === "L") {
      populated = populated || activeSpectra.some(s => s.startsWith("$cl"))
    }
    return populated
  }

  const smartLabel = (label, cats) => {
    const _cats = cats ? cats.split("") : null
    let populated = false
    if (label === "All") {
      populated = Boolean(selectors.filter(i => i.owner).length)
    } else if (label !== "Efficiency") {
      if (activeSpectra.length > 0) {
        populated = _cats.some(c => isPopulated(c))
      }
    }
    return (
      <span className={`tab-header ${populated ? " populated" : ""}`}>
        <span className={classes.bigShow}>
          {label}
          {' '}
          {populated ? " ✶" : ""}
        </span>
        <span className={classes.bigHide}>
          {categoryIcon(_cats && _cats[_cats.length - 1], "", {
            style: { position: "relative", left: 0, height: "1.1rem" },
          })}
        </span>
      </span>
    )
  }

  selectors.sort(selectorSorter)
  const options = useMemo(() => Object.values(ownerInfo), [ownerInfo])

  return (
    <div className="tab-wrapper">
      <Tabs
        value={tab}
        onChange={handleTabChange}
        indicatorColor="primary"
        textColor="primary"
        centered={!ISTOUCH}
        variant={ISTOUCH ? "scrollable" : "standard"}
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
          label={smartLabel("Fluorophores", "DP")}
        />
        <Tab
          className={classes.tabLabel}
          tabIndex={-1}
          label={smartLabel("Filters", "F")}
        />
        <Tab
          tabIndex={-1}
          className={classes.tabLabel}
          label={smartLabel("Light Sources", "L")}
        />
        <Tab
          tabIndex={-1}
          className={classes.tabLabel}
          label={smartLabel("Detectors", "C")}
        />
        <Tab
          tabIndex={-1}
          className={classes.tabLabel}
          label={smartLabel("Efficiency", "%")}
        />
      </Tabs>

      {Object.keys(ownerInfo).length === 0 ? (
        <div className="sweet-loading">
          <RingLoader
            css={override}
            sizeUnit="px"
            size={100}
            color="#ccc"
            loading
          />
        </div>
      ) : (
        <div>
          {tab === 0 && (
            <div>
              <SpectrumSelectorGroup
                selectors={selectors}
                options={options}
                showCategoryIcon
                ownerInfo={ownerInfo}
              />
            </div>
          )}
          {tab === 1 && (
            <div>
              <Typography variant="h6" className={classes.categoryHeader}>
                Fluorescent Proteins
              </Typography>

              <SpectrumSelectorGroup
                selectors={selectors}
                options={options}
                showCategoryIcon
                ownerInfo={ownerInfo}
                category="P"
                hint="protein"
              />

              <Typography variant="h6" className={classes.categoryHeader}>
                Dyes
              </Typography>
              <SpectrumSelectorGroup
                selectors={selectors}
                options={options}
                showCategoryIcon
                ownerInfo={ownerInfo}
                category="D"
                hint="dye"
              />
            </div>
          )}
          {tab === 2 && (
            <div>
              <SpectrumSelectorGroup
                selectors={selectors}
                options={options}
                showCategoryIcon
                ownerInfo={ownerInfo}
                category="F"
                hint="filter"
              />
              <CustomFilterGroup activeSpectra={activeSpectra} />
            </div>
          )}
          {tab === 3 && (
            <div>
              <SpectrumSelectorGroup
                selectors={selectors}
                options={options}
                showCategoryIcon
                ownerInfo={ownerInfo}
                category="L"
                hint="light source"
              />
              <CustomLaserGroup activeSpectra={activeSpectra} />
            </div>
          )}
          {tab === 4 && (
            <div>
              <SpectrumSelectorGroup
                selectors={selectors}
                options={options}
                showCategoryIcon
                ownerInfo={ownerInfo}
                category="C"
                hint="detector"
              />
            </div>
          )}
          {tab === 5 && <EfficiencyTable spectraInfo={spectraInfo} />}
        </div>
      )}
    </div>
  )
})

export default OwnersContainer
