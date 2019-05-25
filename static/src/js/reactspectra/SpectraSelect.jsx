import React, { useState, useContext, useEffect, useRef } from "react";
import Select from "react-select";
import Box from "@material-ui/core/Box";
import ToggleButton from "@material-ui/lab/ToggleButton";
import ToggleButtonGroup from "@material-ui/lab/ToggleButtonGroup";
import IconButton from "@material-ui/core/IconButton";
import Button from "@material-ui/core/Button";
import DeleteIcon from "@material-ui/icons/Delete";
import AddIcon from "@material-ui/icons/Add";
import Tabs from "@material-ui/core/Tabs";
import Tab from "@material-ui/core/Tab";
import { makeStyles } from "@material-ui/core/styles";
import { ID } from "./util";
import { updateSpectra } from "./reducer";
import { AppContext } from "./Store";

const useStyles = makeStyles(theme => ({
  button: {
    "&:hover": {
      // you want this to be the same as the backgroundColor above
      backgroundColor: "#f8f8f8"
    }
  },
  toggleButton: { height: "38px" },
  toggleButtonGroup: { marginLeft: "5px" },
  selector: {
    paddingTop: "6px"
  },
  tabNav: {
    marginBottom: "18px"
  }
}));

const customFilterOption = (option, rawInput) => {
  const words = rawInput.split(" ");
  return words.reduce(
    (acc, cur) => acc && option.label.toLowerCase().includes(cur.toLowerCase()),
    true
  );
};

const SpectrumSelectGroup = ({ category, hint, formState }) => {
  const [selectors, setSelectors] = useState(formState);
  const { state } = useContext(AppContext);
  const classes = useStyles();

  const group = state.owners;
  let options = [];
  if (group) {
    options = Object.keys(group)
      .filter(slug => group[slug].category === category)
      .map(slug => ({
        value: slug,
        label: group[slug].label,
        spectra: group[slug].spectra
      }));
  }

  const addRow = () => setSelectors([...selectors, { id: ID() }]);

  const removeSelector = selID => {
    setSelectors(selectors.filter(({ id }) => id !== selID));
  };

  return (
    <>
      {options.length > 0 ? (
        selectors.map(selector => (
          <div style={{ width: "100%" }} key={selector.id}>
            <Box display="flex">
              <Box flexGrow={1} className={classes.selector}>
                <SpectrumSelector
                  id={selector.id}
                  options={options}
                  defaultValue={selector.value}
                  defaultSubtypes={selector.subtypes}
                />
              </Box>
              {selectors.length > 1 ? (
                <Box>
                  <IconButton
                    aria-label="Delete"
                    color="secondary"
                    className={classes.button}
                    onClick={() => removeSelector(selector.id)}
                  >
                    <DeleteIcon />
                  </IconButton>
                </Box>
              ) : (
                <></>
              )}
            </Box>
          </div>
        ))
      ) : (
        <Select placeholder="Loading..." isLoading />
      )}
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

const formstate = {
  D: {
    owner: ["EX", "EM"]
  }
};

const seriesToFormStates = series => {
  return series.reduce((prev, cur) => {
    const next = { ...prev };
    const { category, owner, subtype } = cur;
    if (!(category in next)) next[category] = [];
    const idx = next[category].findIndex(x => x.value === owner.slug);
    if (idx > -1) {
      next[category][idx].subtypes.push(subtype);
    } else {
      next[category].push({
        id: ID(),
        value: owner.slug,
        label: owner.name,
        subtypes: [subtype]
      });
    }
    return next;
  }, {});
};

const SpectraSelectForm = ({ initial, initialTab }) => {
  const classes = useStyles();
  const [formStates, setFormStates] = useState(seriesToFormStates(initial));
  const [currentTab, setCurrentTab] = useState(initialTab);

  function handleChange(event, newValue) {
    setCurrentTab(newValue);
  }

  useEffect(() => {
    formStates.current = seriesToFormStates(initial);
  }, []); // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <div>
      <Tabs
        value={currentTab}
        onChange={handleChange}
        indicatorColor="primary"
        textColor="primary"
        centered
        className={classes.tabNav}
      >
        <Tab label="Fluorophores" />
        <Tab label="Filters" />
        <Tab label="Light Sources" />
        <Tab label="Cameras" />
        <Tab label="Options" />
      </Tabs>
      <TabContainer index={currentTab}>
        <div>
          <SpectrumSelectGroup
            category="P"
            formState={formStates.P}
            hint="Protein"
          />
          <SpectrumSelectGroup
            category="D"
            formState={formStates.D}
            hint="Dye"
          />
        </div>
        <SpectrumSelectGroup
          category="F"
          formState={formStates.F}
          hint="Filter"
        />
        <SpectrumSelectGroup
          category="L"
          formState={formStates.L}
          hint="Light Source"
        />
        <SpectrumSelectGroup
          category="C"
          formState={formStates.Cchild}
          hint="Camera"
        />
        <div>
          <h4>options</h4>
        </div>
      </TabContainer>
    </div>
  );
};

const TabContainer = ({ index, children }) =>
  children &&
  children.map((child, idx) => (
    <div
      key={idx} // eslint-disable-line react/no-array-index-key
      style={{
        display: index === idx ? "block" : "none"
      }}
    >
      {React.cloneElement(child)}
    </div>
  ));

function dyeCompare(a, b) {
  const TYPE_ORDER = ["AB", "A_2P", "EX", "EM"];
  const upperA = a.subtype.toUpperCase();
  const upperB = b.subtype.toUpperCase();

  if (TYPE_ORDER.indexOf(upperA) > TYPE_ORDER.indexOf(upperB)) return 1;
  return -1;
}

const SpectrumSelector = ({ options, defaultValue, defaultSubtypes }) => {
  const [value, setValue] = useState(defaultValue);
  const [subtypes, setSubtypes] = useState(defaultSubtypes);
  const [availableSubtypes, setAvailableSubtypes] = useState([]);
  const { state, dispatch } = useContext(AppContext);

  const { owners, spectra } = state;
  // when the selector changes
  const handleChange = e => {
    // if it's the same as the previous value do nothing
    const newValue = e && e.value;
    if (newValue === value) return;
    setSubtypes([])
    setValue(newValue);
    // otherwise, set the state (so we can know the previous value on change)
    // then update the AppContext
    updateSpectra({
      add: newValue && owners[newValue].spectra,
      remove: value && owners[value].spectra,
      dispatch
    });
  };

  // when the subtype radio buttons are changed
  const handleType = (event, newTypes) => {
    const removed = [...subtypes].filter(x => !newTypes.includes(x));
    const added = [...newTypes].filter(x => !subtypes.includes(x));
    if (removed.length) {
      const remove = removed.map(
        r => availableSubtypes.find(item => item.subtype === r).id
      );
      updateSpectra({ remove, dispatch });
    }
    if (added.length) {
      const add = added.map(
        r => availableSubtypes.find(item => item.subtype === r).id
      );
      updateSpectra({ add, dispatch });
    }
    setSubtypes(newTypes);
  };

  // clean up on unmount
  useEffect(() => {
    return () => {
      updateSpectra({
        add: null,
        remove: value && owners[value].spectra,
        dispatch
      });
    };
  }, [value, dispatch, owners]);

  // get available subtypes on mount
  useEffect(() => {
    if (value && owners[value].spectra.length > 1) {
      const spectraChoices = owners[value].spectra;
      const available = spectraChoices.map(id => ({
        id,
        subtype: spectra[id].subtype
      }));
      available.sort(dyeCompare);
      setAvailableSubtypes(available);
      if (subtypes === null || subtypes.length === 0) {
        setSubtypes(available.map(({ subtype }) => subtype));
      }
    } else {
      setAvailableSubtypes([]);
    }
  }, [value]); // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <>
      <Box display="flex">
        <Box flexGrow={1}>
          <Select
            isClearable
            defaultValue={
              value ? { value, label: owners[value].label } : undefined
            }
            placeholder="Type to search..."
            filterOption={customFilterOption}
            onChange={handleChange}
            options={options}
          />
        </Box>
        {availableSubtypes.length > 1 && (
          <SubtypeSelector
            subtypes={availableSubtypes}
            onChange={handleType}
            value={subtypes}
          />
        )}
      </Box>
    </>
  );
};

SpectrumSelector.defaultProps = {
  defaultSubtypes: []
};

const SubtypeSelector = ({ subtypes, onChange, value }) => {
  const classes = useStyles();

  const available = subtypes.map(({ subtype }) => subtype);
  return (
    <Box>
      <ToggleButtonGroup
        value={value}
        size="small"
        onChange={onChange}
        className={classes.toggleButtonGroup}
      >
        {available.map(st => (
          <ToggleButton
            key={st}
            value={st}
            color="primary"
            className={classes.toggleButton}
          >
            {st.replace(/^A_/g, "")}
          </ToggleButton>
        ))}
      </ToggleButtonGroup>
    </Box>
  );
};

export default SpectraSelectForm;
