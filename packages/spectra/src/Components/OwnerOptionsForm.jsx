import React from "react"
import gql from "graphql-tag"
import { useQuery, useMutation } from "@apollo/react-hooks"
import ToggleButton from "@material-ui/lab/ToggleButton"
import ToggleButtonGroup from "@material-ui/lab/ToggleButtonGroup"
import { Box } from "@material-ui/core"
import { GET_OWNER_OPTIONS } from "../client/queries"

const MUTATE_OWNER_OPTIONS = gql`
  mutation SetExcludeSubtypes($subtypes: [String]) {
    setExcludeSubtypes(excludeSubtypes: $subtypes) @client
  }
`

const OwnerOptionsForm = () => {
  const [setSubtypes] = useMutation(MUTATE_OWNER_OPTIONS)
  const {
    data: { excludeSubtypes },
  } = useQuery(GET_OWNER_OPTIONS)

  const handleChange = (e, subtypes) => {
    setSubtypes({ variables: { subtypes } })
  }
  return (
    <Box style={{ marginTop: 8 }}>
      <span style={{ marginRight: 8 }}>Exclude subtypes:</span>
      <ToggleButtonGroup
        value={excludeSubtypes}
        size="small"
        onChange={handleChange}
      >
        {["AB", "EX", "EM", "2P"].map(st => (
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
        spectrum types selected here will NOT be added by default when adding a
        new fluorophore
      </p>
    </Box>
  )
}

export default OwnerOptionsForm
