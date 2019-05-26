import React, { useState, useContext } from "react";
import Tabs from "@material-ui/core/Tabs";
import Tab from "@material-ui/core/Tab";
import SpectrumSelectGroup from "./SpectrumSelectGroup";
import { AppContext } from "./Store";

const SpectraSelectForm = ({ formStates, initialTab }) => {
  const [currentTab, setCurrentTab] = useState(initialTab);
  const { state } = useContext(AppContext);
  const { ownerCategories } = state;

  function handleChange(event, newValue) {
    setCurrentTab(newValue);
  }

  return (
    <div>
      <Tabs
        value={currentTab}
        onChange={handleChange}
        indicatorColor="primary"
        textColor="primary"
        centered
        style={{ marginBottom: "18px" }}
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
            options={ownerCategories.P}
            formState={formStates.P}
            hint="Protein"
          />
          <SpectrumSelectGroup
            options={ownerCategories.D}
            formState={formStates.D}
            hint="Dye"
          />
        </div>
        <SpectrumSelectGroup
          options={ownerCategories.F}
          formState={formStates.F}
          hint="Filter"
        />
        <SpectrumSelectGroup
          options={ownerCategories.L}
          formState={formStates.L}
          hint="Light Source"
        />
        <SpectrumSelectGroup
          options={ownerCategories.C}
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

export default SpectraSelectForm;
