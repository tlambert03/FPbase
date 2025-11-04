import React from "react"
import { List } from "react-window"

const ITEM_HEIGHT = 35
const MAX_WINDOW_HEIGHT = 300

/**
 * Virtualized MenuList component for react-select using react-window
 * Only renders visible options for better performance with large lists
 */
const WindowedMenuList = ({ children, maxHeight, _getValue }) => {
  const childrenArray = React.Children.toArray(children)

  // Handle empty children
  if (!childrenArray.length) {
    return <div style={{ padding: "8px" }}>{children}</div>
  }

  const height = Math.min(maxHeight || MAX_WINDOW_HEIGHT, childrenArray.length * ITEM_HEIGHT)

  // Row component for react-window 2.x
  const Row = ({ index, style }) => {
    return <div style={style}>{childrenArray[index]}</div>
  }

  return (
    <List
      rowComponent={Row}
      rowCount={childrenArray.length}
      rowHeight={ITEM_HEIGHT}
      rowProps={{}}
      style={{ height, width: "100%" }}
    />
  )
}

export default WindowedMenuList
