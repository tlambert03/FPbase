import DeleteIcon from "@mui/icons-material/Delete"
import Box from "@mui/material/Box"
import IconButton from "@mui/material/IconButton"
import Typography from "@mui/material/Typography"
import { makeStyles } from "@mui/styles"
import React, { useCallback, useMemo } from "react"
import { useSpectraStore } from "../store/spectraStore"
import { categoryIcon } from "./FaIcon"
import SpectrumSelector from "./SpectrumSelector"

export const useStyles = makeStyles((theme) => ({
  deleteButton: {
    padding: "6px 6px",
    marginRight: 2,
    marginLeft: 2,
    [theme.breakpoints.down("xs")]: {
      display: "none",
    },
    "&:hover": {
      // you want this to be the same as the backgroundColor above
      backgroundColor: "#fff",
      color: "#C84064",
    },
  },
  addButton: {
    marginTop: 6,
    marginRight: 6,
  },
  categoryHeader: {
    textTransform: "uppercase",
    fontSize: "small",
    color: "#3F51B5",
    marginTop: ".35rem",
    marginBottom: "0.4rem",
  },
}))

const SpectrumSelectorGroup = React.memo(function SpectrumSelectorGroup({
  selectors,
  options,
  category = "",
  showCategoryIcon,
  _hint = "item",
  ownerInfo,
  // Optional map of extras to render after each category group in the "All" view
  // e.g. { F: <CustomFilterGroup ... />, L: <CustomLaserGroup ... /> }
  categoryExtras = {},
  // Active spectra to check if custom items exist
  activeSpectra = [],
}) {
  const classes = useStyles()
  const allOwners = useMemo(() => selectors.map(({ owner }) => owner), [selectors])

  // Convert to Set for O(1) lookup in SpectrumSelector (performance optimization)
  const disabledOwners = useMemo(() => new Set(allOwners), [allOwners])

  let mySelectors
  if (category) {
    // For specific categories, include selectors that match the category OR are empty selectors
    mySelectors = selectors.filter((sel) => sel.category === category || !sel.owner)
  } else {
    mySelectors = selectors.filter((sel) => sel.owner || !sel.category)
  }

  let lastCategory = ""
  const categoryNames = {
    P: "Fluorescent Proteins",
    D: "Dyes",
    F: "Filters",
    L: "Light Sources",
    C: "Detectors",
  }

  const categoryOptions = useMemo(
    () => options.filter((opt) => (category ? opt.category === category : true)),
    [category, options]
  )

  const updateActiveSpectra = useSpectraStore((state) => state.updateActiveSpectra)

  const removeRow = useCallback(
    (selector) => {
      // Just remove spectra - selectors will be derived automatically
      if (ownerInfo[selector.owner] && ownerInfo[selector.owner].spectra) {
        const spectraToRemove = ownerInfo[selector.owner].spectra.map(({ id }) => id)
        updateActiveSpectra(undefined, spectraToRemove)
      }
    },
    [ownerInfo, updateActiveSpectra]
  )

  // Track which categories have been rendered to ensure extras are shown
  const renderedCategories = useMemo(() => {
    const cats = new Set(mySelectors.filter((s) => s.owner).map((s) => s.category))
    return cats
  }, [mySelectors])

  // Determine which categories have extras but no selectors
  // Only show if the extras will actually render content
  const categoriesWithExtrasOnly = useMemo(() => {
    if (!category && Object.keys(categoryExtras).length > 0) {
      return Object.keys(categoryExtras).filter((cat) => {
        // Only include category if it has no database selectors
        if (renderedCategories.has(cat)) return false

        // Check if there are custom items for this category
        const prefix = cat === "F" ? "$cf" : cat === "L" ? "$cl" : null
        if (!prefix) return false

        // Only show header if there are active custom items with this prefix
        return activeSpectra.some((id) => id.startsWith(prefix))
      })
    }
    return []
  }, [category, categoryExtras, renderedCategories, activeSpectra])

  return (
    <>
      {mySelectors.map((selector, idx) => {
        const showCategoryHeader = !category && selector.category !== lastCategory
        if (showCategoryHeader) {
          lastCategory = selector.category
        }
        const nextCategory = mySelectors[idx + 1]?.category
        const isLastInCategory = nextCategory !== selector.category
        return (
          <div style={{ width: "100%", margin: "4px 0" }} key={selector.id}>
            {showCategoryHeader && (
              <Typography variant="h6" className={classes.categoryHeader}>
                {categoryNames[selector.category]}
              </Typography>
            )}
            <Box display="flex" alignItems="center" className={classes.root}>
              {categoryIcon(selector.category, "rgba(0,0,50,0.4)", {
                style: {
                  position: "relative",
                  left: category === "L" ? 4 : 2,
                  height: "1.3rem",
                  marginRight: 10,
                },
              })}
              <Box flexGrow={1}>
                <SpectrumSelector
                  key={selector.id}
                  // this line restricts the options to similar categories
                  options={categoryOptions}
                  disabledOwners={disabledOwners}
                  showCategoryIcon={showCategoryIcon}
                  selector={selector}
                  ownerInfo={ownerInfo}
                />
              </Box>
              {selector.owner ? (
                <Box>
                  <IconButton
                    aria-label="Delete"
                    color="secondary"
                    className={classes.deleteButton}
                    tabIndex={-1}
                    onClick={() => removeRow(selector)}
                  >
                    <DeleteIcon />
                  </IconButton>
                </Box>
              ) : null}
            </Box>
            {/* If this is the last selector in a category group, render any provided extras */}
            {isLastInCategory && categoryExtras[selector.category] ? (
              <div style={{ marginTop: 6 }}>{categoryExtras[selector.category]}</div>
            ) : null}
          </div>
        )
      })}
      {/* Render extras for categories that have no database items */}
      {categoriesWithExtrasOnly.map((cat) => (
        <div key={`extras-${cat}`}>
          <Typography variant="h6" className={classes.categoryHeader}>
            {categoryNames[cat]}
          </Typography>
          <div style={{ marginTop: 6 }}>{categoryExtras[cat]}</div>
        </div>
      ))}
    </>
  )
})

export default SpectrumSelectorGroup
