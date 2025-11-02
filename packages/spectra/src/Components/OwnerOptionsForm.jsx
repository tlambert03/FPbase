import { Box } from "@mui/material"
import ToggleButton from "@mui/material/ToggleButton"
import ToggleButtonGroup from "@mui/material/ToggleButtonGroup"
import { useCallback } from "react"
import { useSpectraStore } from "../store/spectraStore"

const OwnerOptionsForm = () => {
  const excludeSubtypes = useSpectraStore((state) => state.excludeSubtypes)
  const setExcludeSubtypes = useSpectraStore((state) => state.setExcludeSubtypes)

  const handleChange = useCallback(
    (_e, subtypes) => {
      setExcludeSubtypes(subtypes)
    },
    [setExcludeSubtypes]
  )
  return (
    <Box style={{ marginTop: 8 }}>
      <span style={{ marginRight: 8 }}>Exclude subtypes:</span>
      <ToggleButtonGroup value={excludeSubtypes} size="small" onChange={handleChange}>
        {["AB", "EX", "EM", "2P"].map((st) => (
          <ToggleButton key={st} value={st} style={{ height: "38px" }}>
            {st}
          </ToggleButton>
        ))}
      </ToggleButtonGroup>
      <p
        style={{
          marginTop: 6,
          fontSize: "small",
          fontStyle: "italic",
          color: "#999",
        }}
      >
        spectrum types selected here will NOT be added by default when adding a new fluorophore
      </p>
    </Box>
  )
}

export default OwnerOptionsForm
