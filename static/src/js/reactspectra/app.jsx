import React, { useEffect, useContext } from "react";
import ReactDOM from "react-dom";
import SpectraSelectForm from "./SpectraSelectForm";
import SpectraViewer from "./SpectraViewer";
import { Store, AppContext } from "./Store";
import LoadingLogo from "./LoadingLogo";
import { ID } from "./util";
// window.onpopstate = event => {};

const stateToSpectraForms = state => {
  const { currentSpectra, spectraInfo } = state;
  return [...new Set(currentSpectra)].reduce((prev, id) => {
    if (!Object.prototype.hasOwnProperty.call(spectraInfo, id)) {
      console.error(`could not get info for spectrum ID ${id}`) // eslint-disable-line
      return prev;
    }
    const { category, owner, label, subtype } = spectraInfo[id];
    if (!Object.prototype.hasOwnProperty.call(prev, category)) {
      // eslint-disable-next-line no-param-reassign
      prev[category] = [];
    }
    const idx = prev[category].findIndex(x => x.value === owner);
    if (idx > -1) {
      prev[category][idx].spectra.push({ id, subtype });
    } else {
      prev[category].push({
        id: ID(),
        value: owner,
        label,
        spectra: [{ id, subtype }]
      });
    }
    return prev;
  }, {});
};

const App = () => {
  const { state } = useContext(AppContext);

  useEffect(() => {
    if (!state.loading) {
      const { currentSpectra } = state;
      let url = window.location.pathname;
      if (currentSpectra.length > 0) {
        url = `?s=${currentSpectra.join(",")}`;
      }
    }
  }, [state]);

  return state && !state.loading ? (
    <div>
      <SpectraViewer />
      <SpectraSelectForm
        formStates={stateToSpectraForms(state)}
        initialTab={state.tab}
      />
    </div>
  ) : (
    <LoadingLogo />
  );
};

const initReactSpectra = elem => {
  ReactDOM.render(
    <Store>
      <App />
    </Store>,
    document.getElementById(elem)
  );
};

export default initReactSpectra;
