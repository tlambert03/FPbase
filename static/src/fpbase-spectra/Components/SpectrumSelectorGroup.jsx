import React from "react"
import Box from "@material-ui/core/Box"
import IconButton from "@material-ui/core/IconButton"
import DeleteIcon from "@material-ui/icons/Delete"
import SpectrumSelector from "./SpectrumSelector"
import { makeStyles } from "@material-ui/core/styles"
import { categoryIcon } from "../Components/FaIcon"
import Typography from "@material-ui/core/Typography"

export const useStyles = makeStyles(theme => ({
  root: {
    [theme.breakpoints.down("sm")]: {
      height: 42
    }
  },
  deleteButton: {
    top: -6,
    [theme.breakpoints.down("sm")]: {
      top: -8
    },
    [theme.breakpoints.down("xs")]: {
      display: "none"
    },
    "&:hover": {
      // you want this to be the same as the backgroundColor above
      backgroundColor: "#fff",
      color: "#C84064"
    }
  },
  addButton: {
    marginTop: 6,
    marginRight: 6
  },
  categoryHeader: {
    textTransform: "uppercase",
    fontSize: "small",
    color: "#3F51B5",
    marginTop: ".2rem",
    marginBottom: "0.4rem"
  }
}))

const SpectrumSelectorGroup = ({
  selectors,
  options,
  category,
  addRow,
  changeOwner,
  removeRow,
  showCategoryIcon,
  hint,
  owners
}) => {
  const classes = useStyles()
  const allOwners = selectors.map(({ owner }) => owner)

  if (category) {
    selectors = selectors.filter(sel => sel.category === category)
  } else {
    selectors = selectors.filter(sel => sel.owner || !sel.category)
  }

  // make sure there is always one empty selector available
  if (selectors.filter(({ owner }) => !owner).length < 1) {
    addRow(category || null)
  }

  let lastCategory = ""
  const categoryNames = {
    P: "Fluorescent Proteins",
    D: "Dyes",
    F: "Filters",
    L: "Light Sources",
    C: "Detectors"
  }

  // disable options that are already claimed by other selectors
  function makeOptions(selector) {
    const otherOwners = allOwners.filter(i => i !== selector.owner)
    return options
      .filter(opt => (category ? opt.category === category : true))
      .map(opt =>
        (otherOwners || []).includes(opt.value)
          ? {
              ...opt,
              label: opt.label + " (already selected)",
              isDisabled: true
            }
          : opt
      )
  }

  return (
    <>
      {selectors.map(selector => (
        <div style={{ width: "100%" }} key={selector.id}>
          {!category &&
            selector.category !== lastCategory &&
            ((lastCategory = selector.category) && (
              <Typography variant="h6" className={classes.categoryHeader}>
                {categoryNames[selector.category]}
              </Typography>
            ))}
          <Box display="flex" className={classes.root}>
            {categoryIcon(selector.category, "rgba(0,0,50,0.4)", {
              style: {
                position: "relative",
                top: 8,
                left: category === "L" ? 4 : 2,
                height: "1.3rem",
                marginRight: 10
              }
            })}
            <Box flexGrow={1}>
              <SpectrumSelector
                key={selector.id}
                // this line restricts the options to similar categories
                options={makeOptions(selector)}
                showCategoryIcon={showCategoryIcon}
                value={selector.owner && owners[selector.owner]}
                onChange={changeOwner(selector.id, category)}
              />
            </Box>
            {(selector.category ? (
              selectors.length > 1
            ) : (
              selectors.filter(({ category }) => !category).length > 1
            )) ? (
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
}

SpectrumSelectorGroup.defaultProps = {
  hint: "item",
  category: ""
}

export default SpectrumSelectorGroup
