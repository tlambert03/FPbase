import React, { memo } from "react"
import Checkbox from "@material-ui/core/Checkbox"
import MUIListItem from "@material-ui/core/ListItem"
import FormControlLabel from "@material-ui/core/FormControlLabel"

const ListCheckbox = memo(function ListCheckbox({
  onCheckItem,
  name,
  value,
  checked,
}) {
  return (
    <MUIListItem>
      <FormControlLabel
        control={
          <Checkbox onChange={onCheckItem} checked={checked} value={value} />
        }
        label={name}
      />
    </MUIListItem>
  )
})

export default ListCheckbox
