module.exports = {
    "env": {
        "es6": true,
        "browser": true,
        "jquery": true,
        "node": true
    },
    "extends": "airbnb",
    "parserOptions": {
       "ecmaVersion": 6
     },
     "rules": {
        "indent": [2, 4],
        "no-unused-vars": 1,
        "object-curly-spacing": "off",
        "require-jsdoc": ["error", {
            "require": {
                "FunctionDeclaration": false,
                "MethodDefinition": false,
                "ClassDeclaration": false,
                "ArrowFunctionExpression": false,
                "FunctionExpression": false
            }
        }]
     }
};
