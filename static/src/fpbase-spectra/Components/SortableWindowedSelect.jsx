import React, { useEffect } from "react"
import { categoryIcon } from "./FaIcon"
import Select, { components } from "react-select"
import { WindowedMenuList } from "react-windowed-select"
import PropTypes from "prop-types"

// const filterOption = ({ label }, query) => {
//   if (!query) return true
//   const words = query.toLowerCase().split(" ")
//   const opt = label.toLowerCase()
//   return words.reduce((acc, cur) => acc && opt.includes(cur), true)
// }

const filterOptions = (query, label) => {
  const words = query.toLowerCase().split(" ")
  const opts = label.toLowerCase()
  return words.reduce((acc, cur) => acc && opts.includes(cur), true)
}

const querySorter = query => {
  const lowerquery = query.trimRight().toLowerCase()
  return function sortOptions(a, b) {
    const alabel = a.label.replace(/^m/, "").toLowerCase()
    const blabel = b.label.replace(/^m/, "").toLowerCase()
    if (alabel.startsWith(lowerquery)) {
      if (blabel.startsWith(lowerquery)) {
        return alabel < blabel ? -1 : 1
      }
      return -1
    }
    if (blabel.startsWith(lowerquery)) return 1
    return alabel < blabel ? -1 : 1
  }
}

const SortableWindowedSelect = ({
  showIcon,
  options,
  components,
  ...otherprops
}) => {
  const [dynamicOptions, setOptions] = React.useState()
  const [inputValue, setInputValue] = React.useState("")

  useEffect(() => {
    setOptions(options)
  }, [options])

  // blur the select element when the escape key is pressed
  const selectRef = React.useRef()
  useEffect(() => {
    const blurme = e => (e.code === "Escape" ? selectRef.current.blur() : null)
    document.addEventListener("keydown", blurme)
    return () => {
      document.removeEventListener("keydown", blurme)
    }
  }, [])

  const handleInputChange = (query, { action }) => {
    setInputValue(query)
    if (action === "input-change") {
      if (query) {
        const newOpts = (options || [])
          .filter(({ label }) => filterOptions(query, label))
          .sort(querySorter(query))
        setOptions(newOpts)
      } else {
        setOptions(options)
      }
    } else if (action === "menu-close") {
      setOptions(options)
    }
  }

  return (
    <Select
      {...otherprops}
      options={dynamicOptions}
      ref={selectRef}
      inputValue={inputValue}
      onInputChange={handleInputChange}
      filterOption={() => true}
      components={{
        ...components,
        ...(showIcon ? { Option: OptionWithIcon } : {}),
        ...{ MenuList: WindowedMenuList }
      }}
    />
  )
}
SortableWindowedSelect.propTypes = {
  showIcon: PropTypes.bool,
  components: PropTypes.object,
  options: PropTypes.array
}

SortableWindowedSelect.defaultProps = {
  showIcon: false
}

const OptionWithIcon = props => {
  const myProps = { ...props }

  myProps.children = (
    <>
      {categoryIcon(props.data.category, "#aaa")}
      {myProps.children}
    </>
  )
  return <components.Option {...myProps} />
}

OptionWithIcon.propTypes = {
  children: PropTypes.node,
  innerProps: PropTypes.object,
  innerRef: PropTypes.oneOfType([PropTypes.func, PropTypes.object]),
  isFocused: PropTypes.bool,
  isSelected: PropTypes.bool
}

export default SortableWindowedSelect
