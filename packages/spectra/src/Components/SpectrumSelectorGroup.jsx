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
  // root: {
  //   [theme.breakpoints.down("sm")]: {
  //     height: 42
  //   }
  // },
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
}) {
  const classes = useStyles()
  const allOwners = useMemo(() => selectors.map(({ owner }) => owner), [selectors])

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

  return (
    <>
      {mySelectors.map((selector) => {
        const showCategoryHeader = !category && selector.category !== lastCategory
        if (showCategoryHeader) {
          lastCategory = selector.category
        }
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
                  allOwners={allOwners}
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
          </div>
        )
      })}
      {/* <Button
        variant="contained"
        color="primary"
        className={classes.addButton}
        onClick={() => addRow(category || null)}
      >
        <AddIcon />
        {`Add ${hint}`}
      </Button> */}
    </>
  )
})

export default SpectrumSelectorGroup
