import React, { useMemo, useCallback } from "react"
import Box from "@material-ui/core/Box"
import IconButton from "@material-ui/core/IconButton"
import DeleteIcon from "@material-ui/icons/Delete"
import { makeStyles } from "@material-ui/core/styles"
import Typography from "@material-ui/core/Typography"
import { useMutation } from "@apollo/react-hooks"
import { categoryIcon } from "./FaIcon"
import SpectrumSelector from "./SpectrumSelector"
import { UPDATE_ACTIVE_SPECTRA, REMOVE_SELECTOR } from "../client/queries"

export const useStyles = makeStyles(theme => ({
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
  category,
  showCategoryIcon,
  hint,
  ownerInfo,
}) {
  const classes = useStyles()
  const allOwners = useMemo(() => selectors.map(({ owner }) => owner), [
    selectors,
  ])

  let mySelectors
  if (category) {
    mySelectors = selectors.filter(sel => sel.category === category)
  } else {
    mySelectors = selectors.filter(sel => sel.owner || !sel.category)
  }

  // make sure there is always one empty selector available
  // if (selectors.filter(({ owner }) => !owner).length < 1) {
  //   addRow(category || null)
  // }

  let lastCategory = ""
  const categoryNames = {
    P: "Fluorescent Proteins",
    D: "Dyes",
    F: "Filters",
    L: "Light Sources",
    C: "Detectors",
  }

  const categoryOptions = useMemo(
    () => options.filter(opt => (category ? opt.category === category : true)),
    [category, options]
  )

  const [removeSelector, { loading: removeLoading }] = useMutation(
    REMOVE_SELECTOR
  )
  const [updateSpectra] = useMutation(UPDATE_ACTIVE_SPECTRA)
  const removeRow = useCallback(
    selector => {
      if (!removeLoading) {
        removeSelector({ variables: { id: selector.id } })
        if (ownerInfo[selector.owner] && ownerInfo[selector.owner].spectra) {
          updateSpectra({
            variables: {
              remove: ownerInfo[selector.owner].spectra.map(({ id }) => id),
            },
          })
        }
      }
    },
    [ownerInfo, removeLoading, removeSelector, updateSpectra]
  )

  return (
    <>
      {mySelectors.map(selector => (
        <div style={{ width: "100%", margin: "4px 0" }} key={selector.id}>
          {!category &&
            selector.category !== lastCategory &&
            ((lastCategory = selector.category) && (
              <Typography variant="h6" className={classes.categoryHeader}>
                {categoryNames[selector.category]}
              </Typography>
            ))}
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
            ) : (
              <></>
            )}
          </Box>
        </div>
      ))}
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

SpectrumSelectorGroup.defaultProps = {
  hint: "item",
  category: "",
}

export default SpectrumSelectorGroup
