import React, { useContext } from "react";
import Box from "@material-ui/core/Box";
import IconButton from "@material-ui/core/IconButton";
import Button from "@material-ui/core/Button";
import DeleteIcon from "@material-ui/icons/Delete";
import AddIcon from "@material-ui/icons/Add";
import { makeStyles } from "@material-ui/core/styles";
import { ID } from "./util";
import SpectrumSelector from "./SpectrumSelector";
import { AppContext } from "./Store";

const useStyles = makeStyles(theme => ({
  button: {
    "&:hover": {
      // you want this to be the same as the backgroundColor above
      backgroundColor: "#f8f8f8"
    }
  }
}));

const SpectrumSelectGroup = ({ category, hint }) => {
  const {
    state: { formState, ownerCategories },
    dispatch
  } = useContext(AppContext);
  const classes = useStyles();

  const options = ownerCategories[category];
  const selectors = formState[category];

  const addRow = () => {
    dispatch({ type: "ADD_FORM_ROW", category });
  };

  const removeRow = id => {
    dispatch({ type: "REMOVE_FORM_ROW", category, id });
  };

  if (!formState[category]) addRow(category);

  if (!(selectors && selectors.length)) {
    return <></>;
  }

  return (
    <>
      {selectors.map(selector => (
        <div style={{ width: "100%" }} key={selector.id}>
          <Box display="flex">
            <Box flexGrow={1} style={{ paddingTop: "6px" }}>
              <SpectrumSelector
                category={category}
                selector={selector}
                options={options}
              />
            </Box>
            {selectors.length > 1 ? (
              <Box>
                <IconButton
                  aria-label="Delete"
                  color="secondary"
                  className={classes.button}
                  onClick={() => removeRow(selector.id)}
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
      <Button
        variant="contained"
        color="primary"
        className="mt-2"
        onClick={addRow}
      >
        {`Add ${hint || "item"}`}
        <AddIcon />
      </Button>
    </>
  );
};

SpectrumSelectGroup.defaultProps = {
  // one id: is required otherwise there will be no rows
  formState: [{ id: ID(), value: null, subtypes: [] }]
};

export default SpectrumSelectGroup;
