import { fetchSpectraList } from "./util";
// Make sure the shape of the default value passed to
// createContext matches the shape that the consumers expect!

const updateSpectra = async ({ add, remove, dispatch }) => {
  if (remove) {
    const toRemove = remove;
    dispatch({ type: "REMOVE_SERIES", toRemove });
  }
  if (add) {
    const toAdd = await fetchSpectraList(add);
    dispatch({ type: "ADD_SERIES", toAdd });
  }
};

function reducer(state, action) {
  switch (action.type) {
    case "REMOVE_SERIES": {
      let newSeries = [...state.series];
      const { toRemove } = action;
      newSeries = newSeries.filter(i => !toRemove.includes(i.id));
      window.series = newSeries;
      return { ...state, series: newSeries };
    }
    case "ADD_SERIES": {
      const series = [...state.series, ...action.toAdd];
      const newState = { ...state, series };
      window.series = newState.series;
      return newState;
    }
    case "UPDATE": {
      const newState = { ...state };
      Object.keys(action).forEach(prop => {
        if (prop !== "type" && prop in newState) {
          newState[prop] = action[prop];
        }
      });
      window.series = newState.series;
      window.owners = newState.owners;
      return newState;
    }
    default:
      throw new Error(`Unrecognized reducer action type: ${action.type}`);
  }
}

export { reducer, updateSpectra };
