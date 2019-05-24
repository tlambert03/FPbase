import React from "react";
// Make sure the shape of the default value passed to
// createContext matches the shape that the consumers expect!
const FilterContext = React.createContext({
  filters: [],
  setFilters: () => {}
});

export { FilterContext };
