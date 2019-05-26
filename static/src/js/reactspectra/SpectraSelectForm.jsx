import React from "react";
import Tabs from "@material-ui/core/Tabs";
import Tab from "@material-ui/core/Tab";
import SpectrumSelectGroup from "./SpectrumSelectGroup";
import { AppContext } from "./Store";

const SpectraSelectForm = () => {
  const {
    state: { tab, formState },
    dispatch
  } = React.useContext(AppContext);

  const handleTabChange = (event, newValue) => {
    if (newValue !== tab) {
      dispatch({ type: "CHANGE_TAB", payload: newValue });
    }
  };

  const isPopulated = category =>
    formState[category] && formState[category].filter(i => i.value).length > 0;

  const smartLabel = (label, cats) => {
    const populated = cats.some(c => isPopulated(c));
    return (
      <span
        style={{
          fontWeight: populated ? "bold" : "normal"
        }}
      >
        {label}
        {populated ? " ✶" : ""}
      </span>
    );
  };

  return (
    <div>
      <Tabs
        value={tab}
        onChange={handleTabChange}
        indicatorColor="primary"
        textColor="primary"
        centered
        style={{ marginBottom: "18px" }}
      >
        <Tab label={smartLabel("Fluorophores", ["D", "P"])} />
        <Tab label={smartLabel("Filters", ["F"])} />
        <Tab label={smartLabel("Light Sources", ["L"])} />
        <Tab label={`Cameras ${isPopulated("C") ? "*" : ""}`} />
        <Tab label="Options" />
      </Tabs>
      <TabContainer index={tab}>
        <div>
          <SpectrumSelectGroup category="P" hint="Protein" />
          <SpectrumSelectGroup category="D" hint="Dye" />
        </div>
        <SpectrumSelectGroup category="F" hint="Filter" />
        <SpectrumSelectGroup category="L" hint="Light Source" />
        <SpectrumSelectGroup category="C" hint="Detector" />
        <div>
          <h4>options</h4>
        </div>
      </TabContainer>
    </div>
  );
};

// here to make sure we always render each tab to maintain the formstate
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

export default SpectraSelectForm;
