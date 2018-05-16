module.exports = {
    "extends": "google",
    "parserOptions": {
       "ecmaVersion": 6
     },
     "rules": {
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
