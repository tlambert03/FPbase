import React, { useEffect } from "react"
import Select, { components } from "react-select"
import { WindowedMenuList } from "react-windowed-select"
import PropTypes from "prop-types"
import { categoryIcon } from "./FaIcon"

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

const OptionWithIcon = props => {
  const myProps = { ...props }

  myProps.children = (
    <>
      {categoryIcon(myProps.data.category, "#aaa")}
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
  isSelected: PropTypes.bool,
}

const emptyFilter = () => true

const SortableWindowedSelect = React.memo(function SortableWindowedSelect({
  showIcon,
  options,
  components,
  ...otherprops
}) {
  const [dynamicOptions, setOptions] = React.useState(options)
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

  const handleInputChange = React.useCallback(
    (query, { action }) => {
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
    },
    [options]
  )

  const memoizedComponents = React.useMemo(
    () => ({
      ...components,
      ...(showIcon ? { Option: OptionWithIcon } : {}),
      ...{ MenuList: WindowedMenuList },
    }),
    [components, showIcon]
  )
  const component = React.useMemo(
    () => (
      <Select
        {...otherprops}
        options={dynamicOptions}
        ref={selectRef}
        inputValue={inputValue}
        onInputChange={handleInputChange}
        filterOption={emptyFilter}
        components={memoizedComponents}
      />
    ),
    [
      dynamicOptions,
      handleInputChange,
      inputValue,
      memoizedComponents,
      otherprops,
    ]
  )
  return component
})

SortableWindowedSelect.propTypes = {
  showIcon: PropTypes.bool,
  components: PropTypes.object,
  options: PropTypes.array,
}

SortableWindowedSelect.defaultProps = {
  showIcon: false,
}

export default SortableWindowedSelect
